use anyhow::{anyhow, bail, Context};
use molar::{
    prelude::*,
    voronoi_cell::{Vector2f, VoronoiCell},
};
use nalgebra::{DVector, SMatrix, SVector};
use rayon::iter::IntoParallelRefMutIterator;
use serde::Deserialize;
use std::{
    collections::{HashMap, HashSet},
    f32::consts::FRAC_PI_2,
    path::{Path, PathBuf},
    sync::Arc,
};

#[cfg(not(test))]
use log::{info, warn}; // Use log crate when building application

#[cfg(test)]
use std::{println as info, println as warn};

mod stats;
pub use stats::*;

mod lipid_molecule;
pub use lipid_molecule::LipidMolecule;

mod lipid_species;
use lipid_species::{LipidSpecies, LipidSpeciesDescr};

mod vmd_visual;
use vmd_visual::VmdVisual;

mod lipid_group;
use lipid_group::LipidGroup;

pub struct Membrane {
    lipids: Vec<LipidMolecule>,
    // Mapping from system resindexes to lipid ids
    resindex_to_id: HashMap<usize, usize>,
    // Lipid groups
    groups: HashMap<String, LipidGroup>,
    // Shared lipid species
    species: Vec<Arc<LipidSpecies>>,
    // Options
    options: MembraneOptions,
    pbox: PeriodicBox,
    // Monolayers. Patches are computer inside the monolayer only
    monolayers: Vec<usize>,
}

#[derive(Deserialize)]
#[serde(default)]
struct MembraneOptions {
    global_normal: Option<Vector3f>,
    order_type: OrderType,
    output_dir: PathBuf,
    max_smooth_iter: usize,
    cutoff: f32,
    lipids: HashMap<String, LipidSpeciesDescr>,
    groups: Vec<String>,
    sel: String,
    // Num of neighbor shells for patch adjustment (default=0 means do not adjust)
    n_shells_patch: usize,
    // Num of neighbor shells for curvature smoothing (default=0 means do not smooth)
    n_shells_smoothing: usize,
}

impl Default for MembraneOptions {
    fn default() -> Self {
        MembraneOptions {
            global_normal: None,
            cutoff: 2.5,
            max_smooth_iter: 1,
            order_type: OrderType::ScdCorr,
            output_dir: ".".into(),
            lipids: Default::default(),
            groups: vec![],
            sel: "all".into(),
            n_shells_patch: 0,
            n_shells_smoothing: 0,
        }
    }
}

impl Membrane {
    pub fn new<K: SerialKind>(source: &Source<K>, optstr: &str) -> anyhow::Result<Self> {
        // Load options
        info!("Processing membrane options...");
        let options: MembraneOptions = toml::from_str(optstr)?;

        // Create working selection
        let source_sel = source.select(&options.sel)?;

        let mut lipids = vec![];
        let mut species = vec![];
        let mut resindex_to_id = HashMap::<usize, usize>::new();

        // Load lipids from provided source
        for (name, descr) in options.lipids.iter() {
            if let Ok(lips) = source_sel.subsel(&descr.whole) {
                // Split into individual lipids (MutableParallel selections)
                let mut lips = lips.split_par_contig(|p| Some(p.atom.resindex))?;
                // Unwrap all lipids in parallel
                lips.par_iter_mut()
                    .try_for_each(|lip| lip.unwrap_simple())?;

                // Now we need subselections, which are not allowed for MutableParallel
                // So convert to ImmutableParallel
                let lips = unsafe { lips.into_immutable() };

                // Use first lipid to create lipid species object
                info!("Creating {} '{}' lipids", lips.len(), name);
                let sp = Arc::new(LipidSpecies::new(name.clone(), descr.clone(), &lips[0])?);
                species.push(Arc::clone(&sp));

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

                    // Compute markers
                    let head_marker = head_sel.center_of_mass()?;
                    let mid_marker = mid_sel.center_of_mass()?;
                    let tail_marker = tail_end_sel.center_of_mass()?;

                    let mut order = Vec::with_capacity(sp.tails.len());
                    for t in &sp.tails {
                        order.push(DVector::from_element(t.bond_orders.len() - 1, 0.0));
                    }

                    let tail_head_vec = tail_marker - head_marker;

                    // Ids are sequencial from zero
                    let id = lipids.len();
                    // Save mapping to real resindexes in order to be able to add lipids
                    // to groups by resindex
                    resindex_to_id.insert(lip.first_atom().resindex, id);

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
                        order,
                        tail_head_vec,
                        id,
                        valid: true,
                        // Rest is set to default
                        patch_ids: vec![],
                        fitted_patch_points: vec![],
                        neib_ids: vec![],
                        voro_vertexes: vec![],
                        normal: Default::default(),
                        mean_curv: 0.0,
                        gaussian_curv: 0.0,
                        princ_dirs: Default::default(),
                        princ_curvs: Default::default(),
                        area: 0.0,
                    });
                }
            } else {
                // No such lipid species found
                warn!("No '{name}' lipids found");
            }
        }

        // If there are any groups in the input toml, create them
        let mut groups: HashMap<String, LipidGroup> = Default::default();
        for gr in &options.groups {
            groups.insert(gr.to_string(), Default::default());
        }

        Ok(Self {
            lipids,
            groups,
            species,
            options,
            pbox: source.require_box()?.clone(),
            monolayers: vec![],
            resindex_to_id,
        })
    }

    pub fn iter_valid_lipids(&self) -> impl Iterator<Item = &LipidMolecule> + Clone {
        self.lipids.iter().filter(|l| l.valid)
    }

    pub fn iter_all_lipids(&self) -> impl Iterator<Item = &LipidMolecule> + Clone {
        self.lipids.iter()
    }

    fn iter_valid_lipids_mut(&mut self) -> impl Iterator<Item = &mut LipidMolecule> {
        self.lipids.iter_mut().filter(|l| l.valid)
    }

    pub fn with_groups(mut self, group_names: &[impl AsRef<str>]) -> Self {
        for gr in group_names {
            self.groups
                .insert(gr.as_ref().to_string(), Default::default());
        }
        self
    }

    pub fn with_global_normal(mut self, global_normal: impl Into<Option<Vector3f>>) -> Self {
        self.options.global_normal = global_normal.into();
        self
    }

    pub fn with_order_type(mut self, order_type: OrderType) -> Self {
        self.options.order_type = order_type;
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

        self.options.output_dir = path.to_path_buf();
        info!("will use output directory '{}'", path.display());
        Ok(self)
    }

    pub fn with_cutoff(mut self, cutoff: f32) -> Self {
        self.options.cutoff = cutoff;
        self
    }

    pub fn with_max_neib_shells(mut self, n: usize) -> Self {
        self.options.n_shells_patch = n;
        self
    }

    pub fn with_max_iter(mut self, max_iter: usize) -> Self {
        self.options.max_smooth_iter = max_iter;
        self
    }

    pub fn reset_groups(&mut self) {
        // Reset lipid indexes in groups
        for gr in self.groups.values_mut() {
            gr.lipid_ids.clear();
            gr.stats.per_species.clear();
        }
    }

    pub fn reset_valid_lipids(&mut self) {
        for l in self.lipids.iter_mut() {
            l.valid = true;
        }
    }

    // pub fn compute_monolayers(&mut self, d: f32) -> anyhow::Result<()> {
    //     let ind: Vec<(usize, usize)> = distance_search_single_pbc(
    //         d,
    //         self.iter_valid_lipids().map(|l| &l.head_marker),
    //         self.iter_valid_lipids().map(|l| l.id),
    //         self.lipids[0].sel.require_box()?,
    //         PBC_FULL,
    //     );

    //     let conn = SearchConnectivity::from_iter(ind.into_iter());

    // }

    pub fn add_ids_to_group(
        &mut self,
        gr_name: impl AsRef<str>,
        ids: &Vec<usize>,
    ) -> anyhow::Result<()> {
        let gr_name = gr_name.as_ref();
        if let Some(gr) = self.groups.get_mut(gr_name) {
            for id in ids {
                // Check if lipid id is in range
                if *id >= self.lipids.len() {
                    bail!(
                        "lipid id {} is out of bounds 0:{}",
                        id,
                        self.lipids.len() - 1
                    );
                }
                // Check if this lipid is valid
                if !self.lipids[*id].valid {
                    continue;
                }
                // Add this lipid id
                gr.lipid_ids.push(*id);
                // Initialize statistics for this lipid species if not yet done
                let sp = &self.lipids[*id].species.name;
                if !gr.stats.per_species.contains_key(sp) {
                    gr.stats.per_species.insert(
                        sp.to_string(),
                        SpeciesStats::new(
                            &self.lipids[*id].species,
                            self.species.iter().map(|sp| sp.name.to_owned()),
                        ),
                    );
                }
            }
        } else {
            bail!("No such group: {gr_name}");
        }

        Ok(())
    }

    pub fn add_resindeces_to_group(
        &mut self,
        gr_name: impl AsRef<str>,
        resindeces: &Vec<usize>,
    ) -> anyhow::Result<()> {
        let ids = resindeces
            .iter()
            .map(|r| {
                self.resindex_to_id
                    .get(r)
                    .cloned()
                    .ok_or_else(|| anyhow!("no lipid with resindex {r}"))
            })
            .collect::<anyhow::Result<Vec<usize>>>()?;
        self.add_ids_to_group(gr_name, &ids)
    }

    pub fn iter_group(
        &self,
        gr_name: impl AsRef<str>,
    ) -> anyhow::Result<impl Iterator<Item = &LipidMolecule>> {
        let gr_name = gr_name.as_ref();
        Ok(self
            .groups
            .get(gr_name)
            .ok_or_else(|| anyhow!("no such group: {gr_name}"))?
            .lipid_ids
            .iter()
            .map(|id| &self.lipids[*id]))
    }

    pub fn iter_group_valid(
        &self,
        gr_name: impl AsRef<str>,
    ) -> anyhow::Result<impl Iterator<Item = &LipidMolecule>> {
        let gr_name = gr_name.as_ref();
        Ok(self
            .groups
            .get(gr_name)
            .ok_or_else(|| anyhow!("no such group: {gr_name}"))?
            .lipid_ids
            .iter()
            .filter_map(|id| {
                if self.lipids[*id].valid {
                    Some(&self.lipids[*id])
                } else {
                    None
                }
            }))
    }

    pub fn iter_group_ids(
        &self,
        gr_name: impl AsRef<str>,
    ) -> anyhow::Result<impl Iterator<Item = usize> + '_> {
        let gr_name = gr_name.as_ref();
        Ok(self
            .groups
            .get(gr_name)
            .ok_or_else(|| anyhow!("no such group: {gr_name}"))?
            .lipid_ids
            .iter()
            .cloned())
    }

    pub fn iter_group_ids_valid(
        &self,
        gr_name: impl AsRef<str>,
    ) -> anyhow::Result<impl Iterator<Item = usize> + '_> {
        let gr_name = gr_name.as_ref();
        Ok(self
            .groups
            .get(gr_name)
            .ok_or_else(|| anyhow!("no such group: {gr_name}"))?
            .lipid_ids
            .iter()            
            .cloned()
            .filter(|i| self.lipids[*i].valid))
    }

    pub fn compute(&mut self) -> anyhow::Result<()> {
        // Compute patches
        self.compute_patches(self.options.cutoff)?;

        // Get initial normals
        self.compute_initial_normals();

        // Do iterative smoothing
        let mut iter = 0;
        loop {
            if self.options.n_shells_patch > 0 && iter == 0 {
                // If asked to do more than one nearest neighbor then
                // the very first iteration is used to only get the 1-st nearest neibours
                self.smooth();
                // Now compute n-th neighbors and set them as patch
                self.patches_from_nth_shell(self.options.n_shells_patch);
            }
            self.smooth();
            iter += 1;
            if iter >= self.options.max_smooth_iter {
                break;
            }
        }

        // Compute order for all valid lipids
        for i in 0..self.lipids.len() {
            if !self.lipids[i].valid {
                continue;
            }
            self.lipids[i].compute_order(
                self.options.order_type.clone(),
                self.options.global_normal.as_ref(),
            );
        }

        // If asked perform curvature smoothing
        self.smooth_curvature(self.options.n_shells_smoothing);

        // Update group stats
        for gr in self.groups.values_mut() {
            gr.frame_update(&self.lipids)?;
        }

        Ok(())
    }

    fn compute_initial_normals(&mut self) {
        // Compute tail-to-head vectors for all lipids
        // This is unwrapped with the same lipid
        for lip in self.iter_valid_lipids_mut() {
            lip.tail_head_vec = (lip.head_marker - lip.tail_marker).normalize();
        }

        // First pass - average of tail-to-head distances in each patch
        for i in 0..self.lipids.len() {
            if self.lipids[i].valid {
                self.lipids[i].normal = self.lipids[i]
                    .patch_ids
                    .iter()
                    .filter_map(|l| {
                        if self.lipids[*l]
                            .tail_head_vec
                            .angle(&self.lipids[i].tail_head_vec)
                            <= FRAC_PI_2
                        {
                            Some(&self.lipids[*l].tail_head_vec)
                        } else {
                            None
                        }
                    })
                    // Also add the central lipid
                    .chain(std::iter::once(&self.lipids[i].tail_head_vec))
                    .sum::<Vector3f>()
                    .normalize();
            }
        }

        // Second pass - average of vectors from the first pass in each patch
        for i in 0..self.lipids.len() {
            if self.lipids[i].valid {
                self.lipids[i].normal = self.lipids[i]
                    .patch_ids
                    .iter()
                    .filter_map(|l| {
                        if self.lipids[*l].normal.angle(&self.lipids[i].normal) <= FRAC_PI_2 {
                            Some(&self.lipids[*l].normal)
                        } else {
                            None
                        }
                    })
                    // Also add the central lipid
                    .chain(std::iter::once(&self.lipids[i].normal))
                    .sum::<Vector3f>()
                    .normalize();
            }
        }

        // Correct normal orientations if needed
        // for i in 0..self.lipids.len() {
        //     if self.lipids[i].valid {
        //         if self.lipids[i].normal.angle(&self.lipids[i].tail_head_vec) > FRAC_PI_2 {
        //             self.lipids[i].normal *= -1.0;
        //         }
        //     }
        // }
    }

    pub fn finalize(&self) -> anyhow::Result<()> {
        // Check if output directory exists and create if needed
        let path = &self.options.output_dir;
        if !path.exists() {
            info!("Creating output directory '{}'...", path.display());
            std::fs::create_dir_all(path)
                .with_context(|| format!("creating output directory '{}'", path.display()))?;
        }

        info!("Writing results to directory '{}'", path.display());

        // Write results for groups
        for (gr_name, gr) in &self.groups {
            info!("\tGroup '{gr_name}'...");
            gr.stats
                .save_group_stats(self.options.output_dir.as_path(), gr_name)?;
            gr.stats
                .save_order_to_file(self.options.output_dir.as_path(), gr_name)?;
        }
        Ok(())
    }

    pub fn set_state(
        &mut self,
        st: impl Into<Holder<State, ImmutableParallel>>,
    ) -> anyhow::Result<()> {
        let st: Holder<State, ImmutableParallel> = st.into();
        // Go over all lipids and set their states
        for lip in &mut self.lipids {
            lip.set_state(st.new_ref())?;
        }
        Ok(())
    }

    pub fn get_state(&self) -> Holder<State, ImmutableParallel> {
        self.lipids.first().unwrap().sel.get_state()
    }

    fn compute_patches(&mut self, d: f32) -> anyhow::Result<()> {
        let ind: Vec<(usize, usize)> = distance_search_single_pbc(
            d,
            self.iter_valid_lipids().map(|l| &l.head_marker),
            self.iter_valid_lipids().map(|l| l.id),
            self.lipids[0].sel.require_box()?,
            PBC_FULL,
        );

        // Clear all patches first!
        for lip in self.lipids.iter_mut() {
            lip.patch_ids.clear();
        }
        // Add ids to patches
        for (i, j) in ind {
            self.lipids[i].patch_ids.push(j);
            self.lipids[j].patch_ids.push(i);
        }
        Ok(())
    }

    // Given pre-computed neib_ids compute n-th neighbour shell
    // and set it as patches
    fn patches_from_nth_shell(&mut self, n_neib: usize) {
        if n_neib < 1 {
            return;
        }

        for i in 0..self.lipids.len() {
            if !self.lipids[i].valid {
                continue;
            }

            let mut neib_list: HashSet<usize> = self.lipids[i].neib_ids.iter().cloned().collect();
            for _ in 2..n_neib {
                let old_neib_list = neib_list.clone();
                for n in old_neib_list {
                    neib_list.extend(self.lipids[n].neib_ids.iter().cloned());
                }
            }

            self.lipids[i].patch_ids = neib_list.into_iter().collect();
        }
    }

    pub fn smooth_curvature(&mut self, n_neib: usize) {
        if n_neib < 1 {
            return;
        }

        let mean_curv: Vec<_> = self.iter_all_lipids().map(|l| l.mean_curv).collect();
        let gauss_curv: Vec<_> = self.iter_all_lipids().map(|l| l.gaussian_curv).collect();

        for i in 0..self.lipids.len() {
            if !self.lipids[i].valid {
                continue;
            }

            let mut neib_list: HashSet<usize> = self.lipids[i].neib_ids.iter().cloned().collect();
            for _ in 2..n_neib {
                let old_neib_list = neib_list.clone();
                for n in old_neib_list {
                    neib_list.extend(self.lipids[n].neib_ids.iter().cloned());
                }
            }

            let m = neib_list.iter().map(|id| mean_curv[*id]).sum::<f32>();
            self.lipids[i].mean_curv = (mean_curv[i] + m) / (neib_list.len() + 1) as f32;
            let g = neib_list.iter().map(|id| gauss_curv[*id]).sum::<f32>();
            self.lipids[i].gaussian_curv = (gauss_curv[i] + g) / (neib_list.len() + 1) as f32;
        }
    }

    pub fn write_vmd_visualization(&self, fname: impl AsRef<Path>) -> anyhow::Result<()> {
        let mut vis = VmdVisual::new();

        for i in 0..self.lipids.len() {
            let lip = &self.lipids[i];
            if !lip.valid {
                continue;
            }
            // Initial lipid marker
            vis.sphere(&self.lipids[i].head_marker, 0.8, "white");
            // tail-head vector
            vis.arrow(&lip.head_marker, &self.lipids[i].tail_head_vec, "yellow");

            // Fitted lipid marker
            vis.sphere(&lip.head_marker, 0.8, "red");
            // fitted normal
            vis.arrow(&lip.head_marker, &lip.normal, "orange");

            // Voronoi cell
            let n = lip.voro_vertexes.len();
            for i in 0..n {
                let p1 = lip.voro_vertexes[i];
                let p2 = lip.voro_vertexes[(i + 1) % n];
                vis.cylinder(&p1, &p2, "green");
            }
            //vis.cylinder(&lip.voro_vertexes[n - 1], &lip.voro_vertexes[0], "green");

            // Fitted patch
            for p in &lip.fitted_patch_points {
                vis.sphere(p, 0.3, "green");
            }
        }

        vis.save_to_file(fname)?;

        Ok(())
    }

    fn smooth(&mut self) {
        // Save current positions of markers to randomly access from the parallel loop
        let saved_markers = self
            .iter_all_lipids()
            .map(|l| l.head_marker)
            .collect::<Vec<_>>();

        self.lipids
            .par_iter_mut()
            .filter(|lip| lip.valid)
            .for_each(|lip| {
                // Get local-to-lab transform
                let to_lab = lip.get_to_lab_transform();
                // Inverse transform
                let to_local = to_lab.try_inverse().unwrap();
                // Local points
                // Local patch could be wrapped over pbc, so we need to unwrap all neighbors.
                // Central point is assumed to be at local zero.
                let p0 = lip.head_marker;
                let local_points = lip
                    .patch_ids
                    .iter()
                    .map(|j| to_local * self.pbox.shortest_vector(&(saved_markers[*j] - p0)))
                    .collect::<Vec<_>>();

                // Compute fitted surface coefs
                let c = get_quad_coefs(&local_points);
                if c.is_none() {
                    return;
                }
                let quad_coefs = c.unwrap();

                // Do Voronoi stuff
                let mut vc = VoronoiCell::new(-10.0, 10.0, -10.0, 10.0);
                for j in 0..local_points.len() {
                    let p = local_points[j];
                    vc.add_point(&Vector2f::new(p.x, p.y), lip.patch_ids[j]);
                }

                // Find direct neighbours
                let mut n_vert = 0;
                lip.neib_ids = vc
                    .iter_vertex()
                    .filter_map(|v| {
                        n_vert += 1;
                        let id = v.get_id();
                        if id >= 0 {
                            Some(id as usize)
                        } else {
                            None
                        }
                    })
                    .collect();

                // Check if node is valid (there are no wall points)
                // If not valid return and don't do extar work
                if lip.neib_ids.len() < n_vert {
                    lip.valid = false;
                    return;
                }

                // Compute curvatures and fitted normal
                lip.compute_curvature_and_normal(&quad_coefs, &to_lab);

                // Project vertexes into the surface and convert to lab space
                lip.voro_vertexes = vc
                    .iter_vertex()
                    .map(|v| {
                        // For now we save only an offset here because
                        // final posiion of the marker is not yet known
                        Pos::from(to_lab * project_to_surf(v.get_pos(), &quad_coefs))
                    })
                    .collect::<Vec<_>>();

                // Compute area.
                // Central point is still in origin since we didn't translate yet.
                // Thus we can just use vertex vectors as sides of triangles in a triangle fan.
                let n = lip.voro_vertexes.len();
                lip.area = 0.0;
                for i in 0..n {
                    lip.area += 0.5
                        * lip.voro_vertexes[i]
                            .coords
                            .cross(&lip.voro_vertexes[(i + 1) % n].coords)
                            .norm();
                }

                // Check if area is too large (if max_area is set for lipid type)
                if lip.species.max_area > 0.0 && lip.area > lip.species.max_area {
                    lip.valid = false;
                    return;
                }

                //Save fitted positions of patch markers
                lip.fitted_patch_points = local_points
                    .iter()
                    .zip(&lip.patch_ids)
                    .map(|(p, id)| {
                        saved_markers[*id]
                            + to_lab * Vector3f::new(0.0, 0.0, z_surf(p.x, p.y, &quad_coefs) - p.z)
                    })
                    .collect();

                // Update fitted marker
                // local z coord of this point is: z_local= a*x^2 + b*y^2 + c*xy + d*x + e*y + f,
                // but x=y=0, so z_local = f = quad_coefs[5]
                // If local height is too large mark point as invalid
                if quad_coefs[5].abs() > 0.5 {
                    lip.valid = false;
                    return;
                }
                lip.head_marker += to_lab * Vector3f::new(0.0, 0.0, quad_coefs[5]);
            }); //nodes

        // Smooth
        // We start from fitted markers themselves
        let mut smooth_n = vec![1.0; self.lipids.len()];
        let mut smooth_p = self
            .lipids
            .iter()
            .map(|l| l.head_marker)
            .collect::<Vec<_>>();
        // Add projected patch points
        for lip in self.iter_valid_lipids() {
            for (id, p) in lip.patch_ids.iter().zip(lip.fitted_patch_points.iter()) {
                smooth_n[*id] += 1.0;
                smooth_p[*id] += 1.0 * p.coords;
            }
        }
        // Compute averages
        for i in 0..self.lipids.len() {
            if self.lipids[i].valid {
                self.lipids[i].head_marker = smooth_p[i] / smooth_n[i];
            }
        }

        // Now compute actual positions of the Voronoi vertices by adding
        // actual position of the marker
        for lip in self.iter_valid_lipids_mut() {
            for v in &mut lip.voro_vertexes {
                *v += lip.head_marker.coords;
            }
        }

        //self.triangulate();
    }
}

pub fn get_quad_coefs<'a>(local_points: &'a Vec<Vector3f>) -> Option<SVector<f32, 6>> {
    //============================
    // Fitting polynomial
    //============================

    // We fit with polynomial fit = a*x^2 + b*y^2 + c*xy + d*x + e*y + f
    // Thus we need a linear system of size 6
    let mut m = SMatrix::<f32, 6, 6>::zeros();
    let mut rhs = SVector::<f32, 6>::zeros(); // Right hand side and result

    let mut powers = SVector::<f32, 6>::zeros();
    powers[5] = 1.0; //free term, the same everywhere
    for loc in local_points {
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
    m.solve_lower_triangular(&rhs)
}

fn project_to_surf(p: Vector2f, coefs: &SVector<f32, 6>) -> Vector3f {
    // local z coord of point is: z_local= a*x^2 + b*y^2 + c*xy + d*x + e*y + f,
    Vector3f::new(p.x, p.y, z_surf(p.x, p.y, coefs))
}

fn z_surf(x: f32, y: f32, coefs: &SVector<f32, 6>) -> f32 {
    // local z coord of point is: z_local= a*x^2 + b*y^2 + c*xy + d*x + e*y + f,
    let z = coefs[0] * x * x
        + coefs[1] * y * y
        + coefs[2] * x * y
        + coefs[3] * x
        + coefs[4] * y
        + coefs[5];
    z
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, io::Read, path::PathBuf};

    use molar::prelude::*;

    use crate::{LipidSpecies, LipidSpeciesDescr, Membrane};

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
            max_area: 0.0,
        };
        let lip_sp = LipidSpecies::new("DOPE".to_owned(), descr, &pope)?;
        println!("{lip_sp:?}");
        Ok(())
    }

    // #[test]
    // fn test_descr_from_itp() -> anyhow::Result<()> {
    //     let top = FileHandler::open("../../tests/POPE.itp")?.read_topology()?;
    //     let n = top.len();
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
        let n = top.len();
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

    #[allow(dead_code)]
    #[derive(Clone, Debug, serde::Deserialize)]
    struct TestInp {
        cutoff: f32,
        lipids: HashMap<String, LipidSpeciesDescr>,
    }

    #[test]
    fn test_toml2() -> anyhow::Result<()> {
        let toml: TestInp = toml::from_str(
            r#"
            cutoff = 2.5
            [lipids.POPC]
            whole = "resname POPE"
            head = "name P N"
            mid = "name C21 C22"
            tails = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]

            [lipids.POPE]
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
        let path = PathBuf::from("tests");
        let src = Source::serial_from_file(path.join("membr.gro"))?;

        // let path = PathBuf::from("/home/semen/work/Projects/Misha/PG_flipping");
        // let src = Source::serial_from_file(path.join("inp_7.gro"))?;

        let z0 = 5.6; //src.select("not resname TIP3 POT CLA /A+/ /B+/")?.center_of_mass()?.z;
        let mut toml = String::new();
        std::fs::File::open("data/inp.toml")?.read_to_string(&mut toml)?;
        let mut memb = Membrane::new(&src, &toml)?;
        //.with_groups(&["upper", "lower"]);
        // .with_global_normal(Vector3f::z())
        // .with_order_type(OrderType::ScdCorr)
        // .with_output_dir("../target/membr_test_results")?
        // .with_cutoff(2.5)
        // .with_max_iter(1);

        let mut upper = vec![];
        let mut lower = vec![];
        for lip in memb.iter_valid_lipids() {
            if lip.head_marker.z > z0 + 1.0 {
                upper.push(lip.id);
            } else if lip.head_marker.z < z0 - 1.0 {
                lower.push(lip.id);
            }
        }

        memb.add_ids_to_group("upper", &upper)?;
        memb.add_ids_to_group("lower", &lower)?;

        let traj = FileHandler::open(path.join("traj_comp.xtc"))?;
        for st in traj.into_iter().take(10) {
            println!(">> {}", st.get_time());
            memb.set_state(st)?;
            memb.compute()?;
        }

        memb.finalize()?;

        Ok(())
    }

    #[test]
    fn test_vmd_vis() -> anyhow::Result<()> {
        let path = PathBuf::from("tests");
        let src = Source::serial_from_file(path.join("membr.gro"))?;
        let z0 = 5.6; //src.select("not resname TIP3 POT CLA")?.center_of_mass()?.z;
        let mut toml = String::new();
        std::fs::File::open("data/inp.toml")?.read_to_string(&mut toml)?;

        let mut memb = Membrane::new(&src, &toml)?
            .with_max_iter(1)
            // .with_global_normal(Vector3f::z())
            // .with_order_type(OrderType::ScdCorr)
            .with_output_dir("../target/membr_test_results")?;

        let mut upper = vec![];
        let mut lower = vec![];
        for lip in memb.iter_valid_lipids() {
            if lip.head_marker.z > z0 + 1.0 {
                upper.push(lip.id);
            } else if lip.head_marker.z < z0 - 1.0 {
                lower.push(lip.id);
            }
            if lip.sel.first_atom().resname == "POGL" {
                println!("{}", lip.head_marker);
            }
        }

        memb.add_ids_to_group("upper", &upper)?;
        memb.add_ids_to_group("lower", &lower)?;

        memb.compute()?;
        memb.finalize()?;
        memb.write_vmd_visualization("../target/vis.tcl")?;

        Ok(())
    }

    #[test]
    fn test_vmd_vis_cg() -> anyhow::Result<()> {
        let path = PathBuf::from("tests");
        let src = Source::serial_from_file(path.join("cg.gro"))?;
        src.select_all()?.set_same_mass(1.0);

        let mut toml = String::new();
        std::fs::File::open(path.join("cg.toml"))?.read_to_string(&mut toml)?;

        let mut memb =
            Membrane::new(&src, &toml)?.with_output_dir("../target/membr_cg_test_results")?;

        let upper = (0..4225).collect();
        let lower = (4225..4225 * 2).collect();

        memb.add_ids_to_group("upper", &upper)?;
        memb.add_ids_to_group("lower", &lower)?;

        memb.compute()?;
        memb.finalize()?;
        memb.write_vmd_visualization("../target/vis_cg.tcl")?;

        for l in memb.iter_valid_lipids_mut() {
            l.sel.set_same_bfactor(l.mean_curv);
        }
        src.save("../target/colored.pdb")?;

        Ok(())
    }
}
