use molar::prelude::*;
use nalgebra::{DVector, SMatrix, SVector};
use std::sync::Arc;

use crate::lipid_species::LipidSpecies;

#[derive(Debug, Default)]
pub struct SingleLipidProperties {
    pub(super) area: f32,
    pub(super) tilt: f32,
    pub(super) order: Vec<DVector<f32>>,

    pub(super) tail_head_vec: Vector3f,
    pub(super) init_normal: Vector3f,
    pub(super) fitted_normal: Vector3f,
    pub(super) fitted_marker: Pos,
    pub(super) mean_curv: f32,
    pub(super) gaussian_curv: f32,
    pub(super) princ_dirs: SMatrix<f32, 3, 2>,
    pub(super) princ_curvs: SVector<f32, 2>,

    pub(super) patch: Vec<usize>,
}

impl SingleLipidProperties {
    pub fn new(species: &LipidSpecies) -> Self {
        let mut order = Vec::with_capacity(species.tails.len());
        for t in &species.tails {
            order.push(DVector::from_element(t.bond_orders.len() - 1, 0.0));
        }
        Self {
            order,
            ..Default::default()
        }
    }

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
        self.fitted_normal = to_lab * Vector3f::new(*d, *e, -1.0).normalize();
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
}

pub struct LipidMolecule {
    pub(super) sel: Sel<MutableSerial>,
    pub(super) species: Arc<LipidSpecies>,
    pub(super) head_sel: Sel<MutableSerial>,
    pub(super) mid_sel: Sel<MutableSerial>,
    pub(super) tail_end_sel: Sel<MutableSerial>,
    pub(super) tail_sels: Vec<Sel<MutableSerial>>,
    pub(super) head_marker: Pos,
    pub(super) mid_marker: Pos,
    pub(super) tail_marker: Pos,

    pub(super) props: SingleLipidProperties,
}

impl LipidMolecule {
    pub fn update_markers(&mut self) -> anyhow::Result<()> {
        self.head_marker = self.head_sel.center_of_mass_pbc()?;
        self.mid_marker = self.mid_sel.center_of_mass_pbc()?;
        self.tail_marker = self.tail_end_sel.center_of_mass_pbc()?;
        Ok(())
    }

    pub fn compute_order(&mut self, order_type: OrderType, normal: &Vector3f) {
        for i in 0..self.tail_sels.len() {
            self.props.order[i] = self.tail_sels[i]
                .lipid_tail_order(
                    order_type.clone(),
                    &vec![normal.clone()],
                    &self.species.tails[i].bond_orders,
                )
                .unwrap();
        }
    }

    pub fn num_tails(&self) -> usize {
        self.tail_sels.len()
    }

    pub fn set_state(&mut self, st: Holder<State, MutableSerial>) -> anyhow::Result<()> {
        self.sel.set_state(st.clone())?;
        self.head_sel.set_state(st.clone())?;
        self.mid_sel.set_state(st.clone())?;
        for t in &mut self.tail_sels {
            t.set_state(st.clone())?;
        }
        self.update_markers()?;
        Ok(())
    }
}
