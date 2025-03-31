use anyhow::{bail, Context};
use molar::prelude::*;
use nalgebra::{zero, Const, DMatrix, Dyn, Matrix, OMatrix, SMatrix, SVector};
use std::fmt::Write;
use std::{
    any,
    collections::HashMap,
    f32::consts::FRAC_PI_2,
    fmt::write,
    path::{Path, PathBuf},
    sync::Arc,
};

#[cfg(not(test))]
use log::{info, warn}; // Use log crate when building application

#[cfg(test)]
use std::{println as info, println as warn};

mod stats;
use stats::{GroupProperties, StatProperties};

mod lipid_molecule;
use lipid_molecule::{LipidMolecule, SingleLipidProperties};

mod lipid_species;
use lipid_species::{LipidSpecies, LipidSpeciesDescr};

pub struct Membrane {
    // source: Source<MutableSerial>,
    //species: Vec<(String, LipidSpecies)>,
    lipids: Vec<LipidMolecule>,
    groups: HashMap<String, LipidGroup>,
    global_normal: Option<Vector3f>,
    order_type: OrderType,
    output_dir: PathBuf,

    // Local patches
    patches: Vec<Vec<usize>>,
}

impl Membrane {
    pub fn new(source: Source<MutableSerial>, defstr: &str) -> anyhow::Result<Self> {
        // Load species descriptions
        let species_descr: HashMap<String, LipidSpeciesDescr> = toml::from_str(defstr)?;

        // Load lipids from provided source
        //let mut species: HashMap<String, Arc<LipidSpecies>> = Default::default();
        let mut lipids = vec![];
        for (name, descr) in species_descr.iter() {
            if let Ok(lips) = source.select(&descr.whole) {
                let lips = lips.split_resindex::<Vec<_>>()?;
                // Use first lipid to create lipid species object
                info!("Creating {} '{}' lipids", lips.len(), name);
                let sp = Arc::new(LipidSpecies::new(name.clone(), descr.clone(), &lips[0])?);

                // Now create individual lipids
                for lip in lips {
                    //let sp = &species[name];
                    let head_sel = lip.subsel(&sp.head_marker_offsets)?;
                    let mid_sel = lip.subsel(&sp.mid_marker_offsets)?;

                    let mut tail_sels = vec![];
                    let mut tail_ends = vec![];
                    for t in sp.tails.iter() {
                        tail_sels.push(lip.subsel(&t.offsets)?);
                        tail_ends.push(*t.offsets.last().unwrap());
                    }
                    let tail_end_sel = lip.subsel(tail_ends)?;

                    // Unwrap the lipid
                    lip.unwrap_simple()?;
                    // Compute markers
                    let head_marker = head_sel.center_of_mass()?;
                    let mid_marker = mid_sel.center_of_mass()?;
                    let tail_marker = tail_end_sel.center_of_mass()?;

                    lipids.push(LipidMolecule {
                        sel: lip, // Selection is moved to the lipid
                        species: Arc::clone(&sp),
                        head_sel,
                        mid_sel,
                        tail_end_sel,
                        tail_sels,
                        head_marker,
                        mid_marker,
                        tail_marker,
                        stats: SingleLipidProperties::new(&sp),
                    });
                }
            } else {
                // No such lipid species found
                warn!("No '{name}' lipids found");
            }
        }
        Ok(Self {
            // source,
            // species,
            lipids,
            groups: Default::default(),
            global_normal: None,
            order_type: OrderType::ScdCorr,
            output_dir: PathBuf::from("membrane_results"),
            patches: vec![],
        })
    }

    pub fn with_global_normal(mut self, global_normal: impl Into<Option<Vector3f>>) -> Self {
        self.global_normal = global_normal.into();
        self
    }

    pub fn with_order_type(mut self, order_type: OrderType) -> Self {
        self.order_type = order_type;
        self
    }

    pub fn with_output_dir(mut self, output_dir: impl AsRef<Path>) -> anyhow::Result<Self> {
        // Check if we can access this directory or create it if needed
        let path = output_dir.as_ref();
        if !path.exists() {
            std::fs::create_dir_all(path).with_context(|| {
                format!("Failed to create output directory '{}'", path.display())
            })?;
        }

        self.output_dir = path.to_path_buf();
        info!("will use output directory '{}'", path.display());
        Ok(self)
    }

    pub fn iter_lipids(&self) -> impl Iterator<Item = (usize, &LipidMolecule)> {
        self.lipids.iter().enumerate()
    }

    pub fn add_group(&mut self, gr_name: impl AsRef<str>, ids: Vec<usize>) -> anyhow::Result<()> {
        let gr_name = gr_name.as_ref();
        for id in &ids {
            if *id >= self.lipids.len() {
                bail!(
                    "lipid id {} is out of bounds 0:{}",
                    id,
                    self.lipids.len() - 1
                );
            }
        }
        if self.groups.contains_key(gr_name) {
            bail!("group '{}' already exists", gr_name);
        } else {
            info!("Creating group '{}' with {} lipids", gr_name, ids.len());
            self.groups
                .insert(gr_name.to_owned(), LipidGroup::new(&self, ids));
        }
        Ok(())
    }

    pub fn compute(&mut self) -> anyhow::Result<()> {
        // Compute patches
        self.compute_patches()?;

        //=========================
        // Compute initial normals
        //=========================

        // Compute tail-to-head vectors for all lipids
        // This is unwrapped with the same lipid
        for lip in &mut self.lipids {
            lip.stats.tail_head_vec = (lip.head_marker - lip.tail_marker).normalize();
        }

        // First pass - average of tail-to-head distances in each patch
        for i in 0..self.lipids.len() {
            self.lipids[i].stats.init_normal = self.patches[i]
                .iter()
                .map(|l| &self.lipids[*l].stats.tail_head_vec)
                .sum::<Vector3f>()
                .normalize();
        }

        // Second pass - average of vectors from the first pass in eac patch
        for i in 0..self.lipids.len() {
            self.lipids[i].stats.init_normal = self.patches[i]
                .iter()
                .map(|l| &self.lipids[*l].stats.init_normal)
                .sum::<Vector3f>()
                .normalize();
        }

        // Correct normal orientations if needed
        for lip in self.lipids.iter_mut() {
            if lip.stats.init_normal.angle(&lip.stats.tail_head_vec) > FRAC_PI_2 {
                lip.stats.init_normal *= -1.0;
            }
        }

        // Initialize fitted_marker and fitted_normal
        for lip in self.lipids.iter_mut() {
            lip.stats.fitted_marker = lip.head_marker;
            lip.stats.fitted_normal = lip.stats.init_normal;
        }

        //=========================
        // Compute Voronoi
        //=========================
        for i in 0..self.lipids.len() {
            self.compute_voronoi(i);
        }

        // Compute order for all lipids
        for lip in &mut self.lipids {
            lip.compute_order(self.order_type.clone(), &self.global_normal.unwrap());
        }
        // Add stats to groups
        for gr in self.groups.values_mut() {
            gr.add_lipid_stats(&self.lipids)?;
        }
        Ok(())
    }

    pub fn smooth(&mut self) {
        // Array of current smoothed points.
        // Initialized from the lipid head groups
        let points = self.lipids.iter().map(|l| l.head_marker).collect::<Vec<_>>();

    }

    pub fn get_to_lab_transform(&mut self, lip_id: usize) -> Matrix3f {
        let mut to_lab = Matrix3f::zeros();
        // Set axes as two perpendicular vectors
        let n = &self.lipids[lip_id].stats.fitted_normal;
        to_lab.set_column(0, &n.cross(&Vector3f::x()));
        to_lab.set_column(1, &n.cross(&to_lab.column(0)));
        to_lab.set_column(2, &-n);
        to_lab
    }

    pub fn get_quad_coefs(&mut self, lip_id: usize, to_lab: &Matrix3f) -> SVector<f32,6> {
        // Inverse transform
        let to_local = to_lab.try_inverse().unwrap();

        //============================
        // Points in local basis
        //============================
        let pbox = self.lipids[0].sel.require_box().unwrap();
        // Local patch could be wrapped over pbc, so we need to unwrap all neighbors.
        // Central point is assumed to be at local zero.
        let p0 = &self.lipids[lip_id].stats.fitted_marker;

        let local_points = self.patches[lip_id]
            .iter()
            .map(|i| to_local*pbox.shortest_vector(&(self.lipids[*i].stats.fitted_marker - p0)))
            .collect::<Vec<_>>();

        //============================
        // Fitting polynomial
        //============================

        // We fit with polynomial fit = a*x^2 + b*y^2 + c*xy + d*x + e*y + f
        // Thus we need a linear system of size 6
        let mut m = nalgebra::SMatrix::<f32, 6, 6>::zeros();
        let mut rhs = nalgebra::SVector::<f32, 6>::zeros(); // Right hand side and result

        let mut powers = nalgebra::SVector::<f32, 6>::zeros();
        powers[5] = 1.0; //free term, the same everywhere
        for loc in local_points.iter()
        {
            // Compute powers
            powers[0] = loc[0] * loc[0]; //xx
            powers[1] = loc[1] * loc[1]; //yy
            powers[2] = loc[0] * loc[1]; //xy
            powers[3] = loc[0]; //x
            powers[4] = loc[1]; //y
                                // Add to the matrix
            m += powers * powers.transpose();
            // rhs
            rhs += powers * loc[2];
        }

        // Now solve and returs coeffs
        m.solve_lower_triangular(&rhs).unwrap()
    }


    pub fn compute_voronoi(&mut self, lip_id: usize) {
        // Get local-to-lab transform
        let to_lab = self.get_to_lab_transform(lip_id);
        // Compute fitted surface coefs (lab-to-local transform computed inside 
        // and not needed outside this function)
        let quad_coefs = self.get_quad_coefs(lip_id, &to_lab);

        //=====================================
        // Compute curvatures and fitted normal
        //=====================================
        compute_curvature_and_normal(&quad_coefs, &to_lab, &mut self.lipids[lip_id].stats);

        // Update fitted marker
        // local z coord of this point is: z_local= a*x^2 + b*y^2 + c*xy + d*x + e*y + f,
        // but x=y=0, so z_local = f = quad_coefs[5]
        self.lipids[lip_id].stats.fitted_marker += to_lab*Vector3f::new(0.0,0.0,quad_coefs[5]);
    }

    pub fn finalize(&self) -> anyhow::Result<()> {
        info!(
            "Writing results to directory '{}'",
            self.output_dir.display()
        );
        // Write results for groups
        for (gr_name, gr) in &self.groups {
            info!("\tGroup '{gr_name}'...");
            gr.stats
                .save_order_to_file(self.output_dir.as_path(), gr_name)?;
        }
        Ok(())
    }

    pub fn set_state(&mut self, st: Holder<State, MutableSerial>) -> anyhow::Result<()> {
        for lip in &mut self.lipids {
            lip.set_state(st.clone())?;
            lip.sel.unwrap_simple()?;
            lip.update_markers()?;
        }
        Ok(())
    }

    fn compute_patches(&mut self) -> anyhow::Result<()> {
        let ind: Vec<(usize, usize)> = distance_search_single_pbc(
            2.0,
            self.lipids.iter().map(|l| &l.mid_marker),
            0..self.lipids.len(),
            self.lipids[0].sel.require_box()?,
            PBC_FULL,
        );
        let conn = LocalConnectivity::from_iter(ind, self.lipids.len());
        self.patches = conn.0;
        Ok(())
    }

    pub fn write_vmd_visualization(&self, fname: impl AsRef<Path>) -> anyhow::Result<()> {
        let mut vis = VmdVisual::new();

        for lip in self.lipids.iter() {
            // Initial lipid marker
            vis.sphere(&lip.head_marker, 0.8, "white");
            // tail-head vector
            vis.arrow(&lip.head_marker, &lip.stats.tail_head_vec, "yellow");
            // initial normal
            vis.arrow(&lip.head_marker, &lip.stats.init_normal, "cyan");
            
            // Fitted lipid marker
            vis.sphere(&lip.stats.fitted_marker, 0.8, "red");
            // fitted normal
            vis.arrow(&lip.stats.fitted_marker, &lip.stats.fitted_normal, "orange");
        }

        vis.save_to_file(fname)?;

        Ok(())
    }
}

//===========================================
#[allow(non_snake_case)]
fn compute_curvature_and_normal(
    coefs: &SVector<f32, 6>,
    to_lab: &SMatrix<f32, 3, 3>,
    out_props: &mut SingleLipidProperties,
) {
    /* Compute the curvatures

    First fundamental form:  I = E du^2 + 2F du dv + G dv^2
    E= r_u dot r_u, F= r_u dot r_v, G= r_v dot r_v

    For us parametric variables (u,v) are just (x,y) in local space.
    Derivatives:
        r_u = {1, 0, 2Ax+Cy+D}
        r_v = {0, 1, 2By+Cx+E}

    In central point x=0, y=0 so:
        r_u={1,0,D}
        r_v={0,1,E}

    Thus: E_ =1+D^2; F_ = D*E; G_ = 1+E^2;

    Second fundamental form: II = L du2 + 2M du dv + N dv2
    L = r_uu dot n, M = r_uv dot n, N = r_vv dot n

    Normal is just  n = {0, 0, 1}
    Derivatives:
        r_uu = {0, 0, 2A}
        r_uv = {0 ,0, C}
        r_vv = {0, 0, 2B}

    Thus: L_ = 2A; M_ = C; N_ = 2B;
    */
    let a = &coefs[0];
    let b = &coefs[1];
    let c = &coefs[2];
    let d = &coefs[3];
    let e = &coefs[4];
    // F is not used;

    let E = 1.0 + d * d;
    let F = d * e;
    let G = 1.0 + e * e;

    let L = 2.0 * a;
    let M = c;
    let N = 2.0 * b;

    //Curvatures:
    out_props.gaussian_curv = (L * N - M * M) / (E * G - F * F);
    out_props.mean_curv = 0.5 * (E * N - 2.0 * F * M + G * L) / (E * G - F * F);

    // Compute normal of the fitted surface at central point
    // dx = 2Ax+Cy+D
    // dy = 2By+Cx+E
    // dz = -1
    // Since we are evaluating at x=y=0:
    // norm = {D,E,1}
    out_props.fitted_normal = to_lab * Vector3f::new(*d, *e, -1.0).normalize();
    // Orientation of the normal could be wrong!
    // Have to be flipped according to lipid orientation later

    /* Principal curvatures
        The principal curvatures k1 and k2 are the eigenvalues
        and the principal directions are eigenvectors
        of the shape operator W:
        W = [I]^-1 * [II]
        W = 1/(EG - F^2) * [E L - F M, E M - F N]
                            [G M - F L, G N - F M]
    */
    let mut W = SMatrix::<f32, 2, 2>::zeros();
    W[(0, 0)] = E * L - F * M;
    W[(0, 1)] = E * M - F * N;
    W[(1, 0)] = G * M - F * L;
    W[(1, 1)] = G * N - F * M;
    W *= 1.0 / (E * G - F * F);
    // W is symmetric despite the equations seems to be not!
    let eig = W.symmetric_eigen();

    out_props.princ_dirs = SMatrix::<f32, 3, 2>::from_columns(&[
        to_lab * Vector3f::new(eig.eigenvectors[(0, 0)], eig.eigenvectors[(1, 0)], 0.0),
        to_lab * Vector3f::new(eig.eigenvectors[(0, 1)], eig.eigenvectors[(1, 1)], 0.0),
    ]);

    out_props.princ_curvs = eig.eigenvalues;
}

//===========================================
struct VmdVisual {
    buf: String,
}

impl VmdVisual {
    fn new() -> Self {
        Self { buf: String::new() }
    }

    fn save_to_file(&self, fname: impl AsRef<Path>) -> anyhow::Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(fname)?;
        writeln!(f, "{}", self.buf)?;
        Ok(())
    }

    fn sphere(&mut self, point: &Pos, radius: f32, color: &str) {
        writeln!(self.buf, "draw color {color}").unwrap();
        writeln!(
            self.buf,
            "draw sphere \"{} {} {}\" radius {radius} resolution 12",
            point.x * 10.0,
            point.y * 10.0,
            point.z * 10.0
        )
        .unwrap();
    }

    fn arrow(&mut self, point: &Pos, dir: &Vector3f, color: &str) {
        use std::fmt::Write;

        const LENGTH: f32 = 5.0;

        let p1 = point * 10.0;
        let p2 = p1 + dir * 0.5 * LENGTH;
        let p3 = p1 + dir * 0.7 * LENGTH;

        writeln!(self.buf, "draw color {color}").unwrap();

        writeln!(
            self.buf,
            "draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2 resolution 12",
            p1.x, p1.y, p1.z, p2.x, p2.y, p2.z
        )
        .unwrap();

        writeln!(
            self.buf,
            "draw cone \"{} {} {}\" \"{} {} {}\" radius 0.4 resolution 12\n",
            p2.x, p2.y, p2.z, p3.x, p3.y, p3.z
        )
        .unwrap();
    }
}

pub struct LipidGroup {
    lipid_ids: Vec<usize>,
    stats: GroupProperties,
}

impl LipidGroup {
    fn new(memb: &Membrane, lipid_ids: Vec<usize>) -> Self {
        let mut stats = GroupProperties::default();
        for id in lipid_ids.iter() {
            let name = memb.lipids[*id].species.name.clone();

            if let Some(el) = stats.per_species.get_mut(&name) {
                el.num_lip += 1;
            } else {
                let sp = &memb.lipids[*id].species;
                stats.per_species.insert(name, StatProperties::new(sp));
            }
        }

        for (sp, stats) in &stats.per_species {
            info!("\t{}: {}", sp, stats.num_lip);
        }
        Self { lipid_ids, stats }
    }

    fn add_lipid_stats(&mut self, all_lipids: &Vec<LipidMolecule>) -> anyhow::Result<()> {
        self.stats
            .add_lipid_stats(self.lipid_ids.iter().map(|id| &all_lipids[*id]))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::{io::Read, path::PathBuf};

    use molar::prelude::*;

    use crate::{lipid_species::PredefinedLipidSpecies, LipidSpecies, LipidSpeciesDescr, Membrane};

    #[test]
    fn test_descr() -> anyhow::Result<()> {
        let src = Source::serial_from_file("tests/membr.gro")?;
        let pope = src.select("resid 60")?;
        let descr = LipidSpeciesDescr {
            whole: "resname POPE".into(),
            head: "name P N".into(),
            mid: "name C21 C22".into(),
            tails: vec![
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218"
                    .into(),
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316".into(),
            ],
        };
        let lip_sp = LipidSpecies::new("DOPE".to_owned(), descr, &pope)?;
        println!("{lip_sp:?}");
        Ok(())
    }

    // #[test]
    // fn test_descr_from_itp() -> anyhow::Result<()> {
    //     let top = FileHandler::open("../../tests/POPE.itp")?.read_topology()?;
    //     let n = top.num_atoms();
    //     let src = Source::new_serial(top.into(), State::new_fake(n).into())?;
    //     let pope = src.select_all()?;
    //     let resname = &pope.first_atom().resname;
    //     let descr = LipidSpeciesDescr::pope();
    //     let lip_sp = LipidSpecies::new(descr, &pope)?;
    //     println!("{lip_sp:?}");
    //     Ok(())
    // }

    #[test]
    fn test_descr_serde() -> anyhow::Result<()> {
        let top = FileHandler::open("tests/POPE.itp")?.read_topology()?;
        let n = top.num_atoms();
        let src = Source::new_serial(top.into(), State::new_fake(n).into())?;

        let descr: LipidSpeciesDescr = toml::from_str(
            r#"
            whole = "resname POPE"
            head = "name P N"
            mid = "name C21 C22"
            tails = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]
        "#,
        )?;

        let lip_sp = LipidSpecies::new("DOPE".to_owned(), descr, &src.select_all()?)?;
        println!("{lip_sp:?}");
        Ok(())
    }

    #[test]
    fn test_toml_descr() -> anyhow::Result<()> {
        let mut toml = String::new();
        std::fs::File::open("data/lipid_species.toml")?.read_to_string(&mut toml)?;
        let descr: PredefinedLipidSpecies = toml::from_str(&toml)?;
        println!("{descr:?}");
        Ok(())
    }

    #[test]
    fn test_toml2() -> anyhow::Result<()> {
        let toml: PredefinedLipidSpecies = toml::from_str(
            r#"
            [POPC]
            whole = "resname POPE"
            head = "name P N"
            mid = "name C21 C22"
            tails = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]

            [POPE]
            whole = "resname POPE"
            head = "name P N"
            mid = "name C21 C22"
            tails = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]
        "#,
        )?;
        println!("{toml:?}");
        Ok(())
    }

    #[test]
    fn test_whole() -> anyhow::Result<()> {
        let path = PathBuf::from("/home/semen/work/Projects/Misha/balanced/not_depleted");
        let src = Source::serial_from_file(path.join("inp_7.gro"))?;
        let z0 = src.select("not resname TIP3 POT CLA")?.center_of_mass()?.z;
        let mut toml = String::new();
        std::fs::File::open("data/lipid_species.toml")?.read_to_string(&mut toml)?;
        let mut memb = Membrane::new(src, &toml)?
            .with_global_normal(Vector3f::z())
            .with_order_type(OrderType::ScdCorr)
            .with_output_dir("../target/membr_test_results")?;

        let mut upper = vec![];
        let mut lower = vec![];
        for (id, lip) in memb.iter_lipids() {
            if lip.head_marker.z > z0 {
                upper.push(id);
            } else {
                lower.push(id);
            }
        }

        memb.add_group("upper", upper)?;
        memb.add_group("lower", lower)?;

        let traj = FileHandler::open(path.join("traj_comp.xtc"))?;
        for st in traj.into_iter().take(10) {
            //println!(">> {}",st.get_time());
            memb.set_state(st.into())?;
            memb.compute()?;
        }

        memb.finalize()?;

        Ok(())
    }

    #[test]
    fn test_vmd_vis() -> anyhow::Result<()> {
        let path = PathBuf::from("/home/semen/work/Projects/Misha/PG_flipping/");
        let src = Source::serial_from_file(path.join("last.gro"))?;
        let z0 = src.select("not resname TIP3 POT CLA")?.center_of_mass()?.z;
        let mut toml = String::new();
        std::fs::File::open("data/lipid_species.toml")?.read_to_string(&mut toml)?;
        let mut memb = Membrane::new(src, &toml)?
            .with_global_normal(Vector3f::z())
            .with_order_type(OrderType::ScdCorr)
            .with_output_dir("../target/membr_test_results")?;

        let mut upper = vec![];
        let mut lower = vec![];
        for (id, lip) in memb.iter_lipids() {
            if lip.head_marker.z > z0 {
                upper.push(id);
            } else {
                lower.push(id);
            }
        }

        memb.add_group("upper", upper)?;
        memb.add_group("lower", lower)?;

        memb.compute()?;
        memb.finalize()?;
        memb.write_vmd_visualization(path.join("vis.tcl"))?;

        Ok(())
    }
}
