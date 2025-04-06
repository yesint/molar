use molar::{
    core::{Matrix3f, PeriodicBox, Pos, Vector3f},
    voronoi_cell::{Vector2f, VoronoiCell},
};
use nalgebra::{SMatrix, SVector};
use rayon::prelude::*;

pub(super) struct Surface {
    pub(super) pbox: PeriodicBox,
    pub(super) nodes: Vec<SurfNode>,
}

#[derive(Default)]
pub(super) struct SurfNode {
    // Computation-related 
    pub(super) valid: bool,
    pub(super) marker: Pos,
    pub(super) patch_ids: Vec<usize>,
    pub(super) fitted_patch_points: Vec<Pos>,
    // Neighbours
    pub(super) neib_ids: Vec<usize>,
    pub(super) voro_vertexes: Vec<Pos>,
    // Curvature
    pub(super) normal: Vector3f,
    pub(super) mean_curv: f32,
    pub(super) gaussian_curv: f32,
    pub(super) princ_dirs: SMatrix<f32, 3, 2>,
    pub(super) princ_curvs: SVector<f32, 2>,
    pub(super) area: f32,
}

impl SurfNode {
    #[allow(non_snake_case)]
    pub(super) fn compute_curvature_and_normal(
        &mut self,
        coefs: &SVector<f32, 6>,
        to_lab: &SMatrix<f32, 3, 3>,
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
        self.gaussian_curv = (L * N - M * M) / (E * G - F * F);
        self.mean_curv = 0.5 * (E * N - 2.0 * F * M + G * L) / (E * G - F * F);

        // Compute normal of the fitted surface at central point
        // dx = 2Ax+Cy+D
        // dy = 2By+Cx+E
        // dz = -1
        // Since we are evaluating at x=y=0:
        // norm = {D,E,1}
        self.normal = to_lab * Vector3f::new(*d, *e, -1.0).normalize();
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

        self.princ_dirs = SMatrix::<f32, 3, 2>::from_columns(&[
            to_lab * Vector3f::new(eig.eigenvectors[(0, 0)], eig.eigenvectors[(1, 0)], 0.0),
            to_lab * Vector3f::new(eig.eigenvectors[(0, 1)], eig.eigenvectors[(1, 1)], 0.0),
        ]);

        self.princ_curvs = eig.eigenvalues;
    }

    pub fn get_to_lab_transform(&self) -> Matrix3f {
        let mut to_lab = Matrix3f::zeros();
        to_lab.set_column(0, &self.normal.cross(&Vector3f::x()));
        to_lab.set_column(1, &self.normal.cross(&to_lab.column(0)));
        to_lab.set_column(2, &-self.normal);
        to_lab
    }    
}

impl Surface {
    pub(super) fn smooth(&mut self) {
        // Save current positions of markers to randomly access from the parallel loop
        let saved_markers = self.nodes.iter().map(|l| l.marker).collect::<Vec<_>>();

        self.nodes
            .par_iter_mut()
            .filter(|node| node.valid)
            .for_each(|node| {
                // Get local-to-lab transform
                let to_lab = node.get_to_lab_transform();
                // Inverse transform
                let to_local = to_lab.try_inverse().unwrap();
                // Local points
                // Local patch could be wrapped over pbc, so we need to unwrap all neighbors.
                // Central point is assumed to be at local zero.
                let p0 = node.marker;
                let local_points = node
                    .patch_ids
                    .iter()
                    .map(|j| to_local * self.pbox.shortest_vector(&(saved_markers[*j] - p0)))
                    .collect::<Vec<_>>();

                // Compute fitted surface coefs
                let quad_coefs = get_quad_coefs(&local_points);

                // Do Voronoi stuff
                let mut vc = VoronoiCell::new(-10.0, 10.0, -10.0, 10.0);
                for j in 0..local_points.len() {
                    let p = local_points[j];
                    vc.add_point(&Vector2f::new(p.x, p.y), node.patch_ids[j]);
                }

                // Find direct neighbours
                let mut n_vert = 0;
                node.neib_ids = vc
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
                if node.neib_ids.len() < n_vert {
                    node.valid = false;
                    return;
                }

                // Compute curvatures and fitted normal
                node.compute_curvature_and_normal(&quad_coefs, &to_lab);

                // Project vertexes into the surface and convert to lab space
                node.voro_vertexes = vc
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
                let n = node.voro_vertexes.len();
                for i in 0..n {
                    node.area += 0.5*node.voro_vertexes[i].coords.cross(&node.voro_vertexes[(i+1)%n].coords).norm();
                }

                //Save fitted positions of patch markers
                node.fitted_patch_points = local_points
                    .iter()
                    .zip(&node.patch_ids)
                    .map(|(p, id)| {
                        saved_markers[*id]
                            + to_lab * Vector3f::new(0.0, 0.0, z_surf(p.x, p.y, &quad_coefs) - p.z)
                    })
                    .collect();

                // Update fitted marker
                // local z coord of this point is: z_local= a*x^2 + b*y^2 + c*xy + d*x + e*y + f,
                // but x=y=0, so z_local = f = quad_coefs[5]
                node.marker += to_lab * Vector3f::new(0.0, 0.0, quad_coefs[5]);
            }); //nodes

        // Smooth
        // We start from fitted markers themselves
        let mut smooth_n = vec![1.0; self.nodes.len()];
        let mut smooth_p = self.nodes.iter().map(|l| l.marker).collect::<Vec<_>>();
        // Add projected patch points
        for lip in &self.nodes {
            for (id, p) in lip.patch_ids.iter().zip(lip.fitted_patch_points.iter()) {
                smooth_n[*id] += 1.0;
                smooth_p[*id] += 1.0 * p.coords;
            }
        }
        // Compute averages
        for i in 0..self.nodes.len() {
            self.nodes[i].marker = smooth_p[i] / smooth_n[i];
        }

        // Now compute actual positions of the Voronoi vertices by adding
        // actual position of the marker
        for node in &mut self.nodes {
            for v in &mut node.voro_vertexes {
                *v += node.marker.coords;
            }
        }

        //self.triangulate();
    }

    // fn triangulate(&self) {
    //     let mut points = HashMap::<[usize;3],Vec<TriangulationPoint>>::new();
    //     let mut np = 0;
    //     for cur in 0..self.nodes.len() {
    //         if !self.nodes[cur].valid {
    //             continue;
    //         }
    //         // Neighbours are in a ccw direction
    //         let n = self.nodes[cur].neib.len();
    //         for i in 0..n {
    //             // For point i+1 the parents are:
    //             let a = cur;
    //             let b = self.nodes[cur].neib[i]; // i -> i+1
    //             let c = self.nodes[cur].neib[(i+1)%n]; // i+1 -> i+2
    //             let mut parents = [a,b,c];
    //             parents.sort();
                
    //             let new_point = TriangulationPoint{
    //                 pos: self.nodes[cur].vertexes[(i+1)%n],
    //                 node_id: cur,
    //                 neib_id: (i+1)%n,
    //                 glob_ind: np,
    //             };

    //             if points.contains_key(&parents) {
    //                 points.get_mut(&parents).unwrap().push(new_point);
    //             } else {
    //                 points.insert(parents, vec![new_point]);
    //             }
                
    //             np += 1;
    //         }
    //     }

    //     for (id,p) in &points {
    //         println!("{:?}: {}",id,p.len());
    //     }


    // }
}

struct TriangulationPoint {
    pos: Pos,
    node_id: usize,
    neib_id: usize,
    glob_ind: usize,
}

pub fn get_quad_coefs<'a>(local_points: &'a Vec<Vector3f>) -> SVector<f32, 6> {
    //============================
    // Fitting polynomial
    //============================

    // We fit with polynomial fit = a*x^2 + b*y^2 + c*xy + d*x + e*y + f
    // Thus we need a linear system of size 6
    let mut m = nalgebra::SMatrix::<f32, 6, 6>::zeros();
    let mut rhs = nalgebra::SVector::<f32, 6>::zeros(); // Right hand side and result

    let mut powers = nalgebra::SVector::<f32, 6>::zeros();
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
    m.solve_lower_triangular(&rhs).unwrap()
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
