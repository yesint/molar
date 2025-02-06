mod holder;
mod kinds;
mod sel;
mod sel_split;
mod source;
mod utils;

pub use holder::*;
pub use kinds::*;
pub use sel::*;
pub use sel_split::*;
pub use source::*;
pub(crate) use utils::*;

use super::{selection_parser::SelectionParserError, BuilderError, PeriodicBoxError};
use crate::io::FileIoError;
use thiserror::Error;

//############################################################
//#  Error enums
//############################################################

#[derive(Error, Debug)]
#[error("topology and state have different sizes ({0},{1})")]
pub struct TopologyStateSizes(usize, usize);

#[derive(Error, Debug)]
pub enum SelectionError {
    #[error("selection parser failed")]
    Parser(#[from] SelectionParserError),

    #[error(transparent)]
    DifferentSizes(#[from] TopologyStateSizes),

    #[error("creating selection from expression {expr_str}")]
    FromExpr {
        expr_str: String,
        source: SelectionIndexError,
    },

    #[error("creating selection from range {first}:{last}")]
    FromRange {
        first: usize,
        last: usize,
        source: SelectionIndexError,
    },

    #[error("invalid local sub-range: {0}:{1}, valid range: 0:{2}")]
    LocalRange(usize, usize, usize),

    #[error("creating selection from vec {first}..{last} of size {size}")]
    FromVec {
        first: usize,
        last: usize,
        size: usize,
        source: SelectionIndexError,
    },

    #[error("index {0} is beyond the allowed range 0:{1}")]
    OutOfBounds(usize, usize),

    #[error("selection index {0}:{1} is outside the source range: 0:{2}")]
    IndexCheck(usize, usize, usize),

    #[error(transparent)]
    FileIo(#[from] FileIoError),

    #[error(transparent)]
    Builder(#[from] BuilderError),

    #[error("can't set incompatible state")]
    SetState,

    #[error("can't set incompatible topology")]
    SetTopology,

    #[error("can't release source: multiple references are active")]
    Release,

    #[error("splitting produced no selections")]
    EmptySplit,

    #[error("selection intersection is empy")]
    EmptyIntersection,
    
    #[error("selection difference is empy")]
    EmptyDifference,
    
    #[error("selection complement is empy")]
    EmptyComplement,

    #[error(transparent)]
    PeriodicBox(#[from] PeriodicBoxError),
}

#[derive(Error, Debug)]
pub enum SelectionIndexError {
    #[error("selection index is empty")]
    IndexEmpty,
    #[error("selection indeces {0}:{1} are out of allowed range 0:{2}")]
    IndexOutOfBounds(usize, usize, usize),
}

//############################################################
//#  Tests
//############################################################

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    pub use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
    use rayon::iter::IndexedParallelIterator;
    use std::sync::LazyLock;

    static TOPST: LazyLock<(Topology, State)> = LazyLock::new(|| {
        let mut h = FileHandler::open("tests/protein.pdb").unwrap();
        h.read().unwrap()
    });

    #[test]
    fn builder_overlap() -> anyhow::Result<()> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = Source::new_serial(top.into(), st.into())?;
        // Create two overlapping selections
        let _sel1 = b.select_iter(0..10)?;
        let _sel2 = b.select_iter(5..15)?;
        Ok(())
    }

    #[test]
    fn builder_par_no_overlap() {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = Source::new_serial(top.into(), st.into()).unwrap();
        // Create two non-overlapping selections.
        let _sel1 = b.select_iter(0..10).unwrap();
        let _sel2 = b.select_iter(11..15).unwrap();
    }

    fn make_sel_all() -> anyhow::Result<Sel<MutableSerial>> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = Source::new_serial(top.into(), st.into())?;
        let sel = b.select_all()?;
        Ok(sel)
    }

    fn make_sel_prot() -> anyhow::Result<Sel<MutableSerial>> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = Source::new_serial(top.into(), st.into())?;
        let sel = b.select_str("not resname TIP3 POT CLA")?;
        Ok(sel)
    }

    #[test]
    fn test_measure() -> anyhow::Result<()> {
        let sel = make_sel_all()?;
        println!("before {}", sel.iter_pos().next().unwrap());

        let (minv, maxv) = sel.min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}", sel.iter_pos().next().unwrap());
        println!("{:?}", sel.min_max());
        Ok(())
    }

    #[test]
    fn test_measure_pbc() -> anyhow::Result<()> {
        let sel = make_sel_all()?;

        let cm = sel.center_of_mass()?;
        println!("{cm}");
        Ok(())
    }

    #[test]
    fn test_translate() -> anyhow::Result<()> {
        let sel = make_sel_all()?;

        println!("before {}", sel.iter_pos().next().unwrap());
        sel.translate(&Vector3f::new(10.0, 10.0, 10.0));
        println!("after {}", sel.iter_pos().next().unwrap());
        Ok(())
    }

    #[test]
    fn test_write_to_file() -> anyhow::Result<()> {
        let sel = make_sel_all()?;

        let mut h = FileHandler::create(concat!(env!("OUT_DIR"), "/f.pdb"))?;
        h.write_topology(&sel)?;
        h.write_state(&sel)?;
        Ok(())
    }

    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let sel = make_sel_prot()?;
        sel.unwrap_connectivity_dim(0.2, PBC_FULL)?;

        let mut h = FileHandler::create(concat!(env!("OUT_DIR"), "/unwrapped.pdb"))?;
        h.write_topology(&sel)?;
        h.write_state(&sel)?;
        Ok(())
    }

    #[test]
    fn eigen_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let sel2 = make_sel_prot()?;

        sel2.rotate(&Vector3f::x_axis(), 80.0_f32.to_radians());

        sel2.save(concat!(env!("OUT_DIR"), "/sel2.pdb"))?;
        sel1.save(concat!(env!("OUT_DIR"), "/sel1_before.pdb"))?;
        println!("Initial RMSD:{}", Sel::rmsd_mw(&sel1, &sel2)?);

        let m = Sel::fit_transform(&sel1, &sel2)?;
        println!("{m}");

        sel1.apply_transform(&m);

        sel1.save(concat!(env!("OUT_DIR"), "/sel1_after.pdb"))?;
        println!("Final RMSD:{}", Sel::rmsd_mw(&sel1, &sel2)?);
        Ok(())
    }

    #[test]
    fn split_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        for res in sel1.split_iter(|p| Some(p.atom.resid))? {
            println!("Res: {}", res.iter_atoms().next().unwrap().resid)
        }
        Ok(())
    }

    #[test]
    fn sasa_test() -> anyhow::Result<()> {
        let sel1 = make_sel_all()?;
        let (a, v) = sel1.sasa();
        println!("Sasa: {a}, Volume: {v}");
        Ok(())
    }

    #[test]
    fn tets_gyration() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let g = sel1.gyration_pbc()?;
        println!("Gyration radius: {g}");
        Ok(())
    }

    #[test]
    fn test_inertia() -> anyhow::Result<()> {
        let sel1 = make_sel_all()?;
        //sel1.rotate(
        //    &UnitVector3::new_normalize(Vector3f::new(1.0,2.0,3.0)),
        //    0.45
        //);
        let (moments, axes) = sel1.inertia_pbc()?;
        println!("Inertia moments: {moments}");
        println!("Inertia axes: {axes:?}");
        Ok(())
        // {0.7308828830718994 0.3332606256008148 -0.5956068634986877}
        // {0.6804488301277161 -0.28815552592277527 0.6737624406814575}
        // {-0.05291106179356575 0.8977214694023132 0.4373747706413269}
    }

    #[test]
    fn test_principal_transform() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let tr = sel1.principal_transform()?;
        println!("Transform: {tr}");

        let (_, axes) = sel1.inertia()?;
        println!("Axes before: {axes}");

        sel1.apply_transform(&tr);

        let (_, axes) = sel1.inertia()?;
        println!("Axes after: {axes}");

        sel1.save(concat!(env!("OUT_DIR"), "/oriented.pdb"))?;

        Ok(())
    }

    #[test]
    fn test_builder_append_from_self() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let n = top.num_atoms();
        let builder = Source::new_builder(top.into(), st.into())?;
        let sel = builder.select_str("resid 550:560")?;
        let added = sel.len();
        builder.append(&sel);
        let all = builder.select_all()?;
        assert_eq!(all.len(), n + added);
        Ok(())
    }

    #[test]
    fn test_builder_remove_from_self() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let n = top.num_atoms();
        let builder = Source::new_builder(top.into(), st.into())?;
        let sel = builder.select_str("resid 550:560")?;
        let removed = sel.len();
        builder.remove(&sel)?;
        let all = builder.select_all()?;
        assert_eq!(all.len(), n - removed);
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_builder_fail_invalid_selection() {
        let (top, st) = FileHandler::open("tests/protein.pdb")
            .unwrap()
            .read()
            .unwrap();
        let builder = Source::new_builder(top.into(), st.into()).unwrap();
        let sel = builder.select_str("resid 809").unwrap(); // last residue
        builder.remove(&sel).unwrap();
        // Trying to call method on invalid selection
        let _cm = sel.center_of_mass();
    }

    // #[test]
    // fn as_parrallel_test() -> anyhow::Result<()> {
    //     let top = TOPST.0.clone();
    //     let st = TOPST.1.clone();
    //     let ser = Source::new_serial(top.into(), st.into())?;
    //     let ser_sel = ser.select_all()?;

    //     let mut res = vec![];

    //     ser_sel.as_parallel(|par_sel| {
    //         let fragments = par_sel.into_iter_fragments_resindex();
    //         res = fragments.par_iter().map(|sel| {
    //             let cm = sel.center_of_mass().unwrap();
    //             sel.translate(&cm.coords);
    //             println!("thread: {}", rayon::current_thread_index().unwrap());
    //             sel.center_of_mass().unwrap()
    //         }).collect::<Vec<_>>();
    //     })?;

    //     println!("cm = {:?}",res);

    //     Ok(())
    // }

    #[test]
    fn as_parrallel_test() -> anyhow::Result<()> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let ser = Source::new_serial(top.into(), st.into())?;
        let ser_sel = ser.select_all()?;

        let mut parts = ser_sel.split_par_contig(|p| Some(p.atom.resindex))?;
        let coms = parts
            .par_iter()
            .map(|sel| sel.center_of_mass())
            .collect::<Result<Vec<_>,_>>()?;

        parts
            .par_iter()
            .enumerate()
            .for_each(|(i, sel)| sel.translate(&(-coms[i].coords)));

        println!("{:?}", coms);

        Ok(())
    }

    #[test]
    fn test_or() -> anyhow::Result<()> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let ser = Source::new_serial(top.into(), st.into())?;
        let sel1 = ser.select_str("name CA")?;
        let sel2 = ser.select_str("name CB")?;
        let sel3 = sel1.union(&sel2);
        assert_eq!(sel3.len(), sel1.len()+sel2.len());
        Ok(())
    }
}
