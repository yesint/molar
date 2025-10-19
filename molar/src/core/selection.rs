mod selection_def;
mod utils;

use std::{cell::{BorrowError, BorrowMutError}, num::ParseIntError};

// pub use holder::*;
// pub use kinds::*;
// pub use sel::*;
// pub use sel_split::*;
// pub use source::*;
pub use molar_powersasa::SasaResults;
pub use selection_def::*;
pub(crate) use utils::*;
mod sel;
pub use sel::*;

use super::{selection_parser::SelectionParserError, BuilderError, PeriodicBoxError};
use crate::io::FileIoError;
use thiserror::Error;

/// Alias to sorted vector
pub type SVec = sorted_vec::SortedSet<usize>;

//############################################################
//#  Error enums
//############################################################

/// Error for different sizes of topology and state
#[derive(Error, Debug)]
#[error("topology and state have different sizes ({0},{1})")]
pub struct TopologyStateSizesError(usize, usize);

/// Error related to creation of selections
#[derive(Error, Debug)]
pub enum SelectionError {
    #[error("can't bind selection")]
    Borrow(#[from] BorrowError),

    #[error("can't bind selection mutably")]
    BorrowMut(#[from] BorrowMutError),

    #[error("selection parser failed")]
    Parser(#[from] SelectionParserError),

    #[error(transparent)]
    DifferentSizes(#[from] TopologyStateSizesError),

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

    #[error("local index {0} is beyond the allowed range 0:{1}")]
    LocalToGlobal(usize, usize),

    #[error("selection index {0}:{1} is outside the source range: 0:{2}")]
    IndexValidation(usize, usize, usize),

    #[error(transparent)]
    FileIo(#[from] FileIoError),

    #[error(transparent)]
    Builder(#[from] BuilderError),

    #[error("can't set incompatible state")]
    IncompatibleState,

    #[error("can't set incompatible topology")]
    IncompatibleTopology,

    #[error("can't release source: multiple references are active")]
    Release,

    #[error("selection from vector slice is empty")]
    EmptySlice,

    #[error("selection range {0}:{1} is empty")]
    EmptyRange(usize, usize),

    #[error("selection '{0}' is empty")]
    EmptyExpr(String),

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

    #[error("no molecules in topology")]
    NoMolecules,

    #[error("gromacs ndx error")]
    Ndx(#[from] NdxError),

    #[error("selection as a definition for subselecting is ambigous")]
    SelDefInSubsel,

    #[error("can't make par_split from overlapping selections")]
    ParSplitOverlap,

    #[error("can't make par_split from different systems")]
    ParSplitDifferentSystems,
}

/// Errors related to accessing selection indexes
#[derive(Error, Debug)]
pub enum SelectionIndexError {
    #[error("selection index is empty")]
    IndexEmpty,
    #[error("selection indeces {0}:{1} are out of allowed range 0:{2}")]
    IndexOutOfBounds(usize, usize, usize),
}

/// Errors related to reading and manipulating Gromacs index files
#[derive(Debug, Error)]
pub enum NdxError {
    #[error("group {0} not found")]
    NoGroup(String),

    #[error("group {0} is empty")]
    EmptyGroup(String),

    #[error("index parse error in group {0}")]
    Parse(String, #[source] ParseIntError),

    #[error("error reading ndx file {0}")]
    NdxIo(std::path::PathBuf, #[source] std::io::Error),

    #[error("malformed ndx file {0}")]
    MalformedNdxFile(std::path::PathBuf),
}

//############################################################
//#  Tests
//############################################################

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use rayon::iter::IndexedParallelIterator;
    pub use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
    use std::sync::LazyLock;

    static TOPST: LazyLock<(Topology, State)> = LazyLock::new(|| {
        let mut h = FileHandler::open("tests/protein.pdb").unwrap();
        h.read().unwrap()
    });

    #[test]
    fn builder_overlap() -> anyhow::Result<()> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = System::new(top, st)?;
        // Create two overlapping selections
        let _sel1 = b.select(0..10)?;
        let _sel2 = b.select(5..15)?;
        Ok(())
    }

    #[test]
    fn builder_par_no_overlap() {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = System::new(top, st).unwrap();
        // Create two non-overlapping selections.
        let _sel1 = b.select(0..10).unwrap();
        let _sel2 = b.select(11..15).unwrap();
    }

    fn make_sel_all() -> anyhow::Result<Sel> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = System::new(top, st)?;
        let sel = b.select("all")?;
        Ok(sel)
    }

    fn make_sel_prot() -> anyhow::Result<Sel> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = System::new(top, st)?;
        let sel = b.select("not resname TIP3 POT CLA")?;
        Ok(sel)
    }

    #[test]
    fn test_measure() -> anyhow::Result<()> {
        let sel = make_sel_all()?;
        let selb = sel.bind()?;
        println!("before {}", selb.iter_pos().next().unwrap());

        let (minv, maxv) = selb.min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}", selb.iter_pos().next().unwrap());
        println!("{:?}", selb.min_max());
        Ok(())
    }

    #[test]
    fn test_measure_pbc() -> anyhow::Result<()> {
        let sel = make_sel_all()?;

        let cm = sel.bind()?.center_of_mass()?;
        println!("{cm}");
        Ok(())
    }

    #[test]
    fn test_translate() -> anyhow::Result<()> {
        let sel = make_sel_all()?;
        let selb = sel.bind()?;

        println!("before {}", selb.iter_pos().next().unwrap());
        sel.bind_mut()?.translate(&Vector3f::new(10.0, 10.0, 10.0));
        println!("after {}", selb.iter_pos().next().unwrap());
        Ok(())
    }

    #[test]
    fn test_write_to_file() -> anyhow::Result<()> {
        let sys = System::from_file("tests/protein.pdb")?;
        let sel = sys.select("name CA")?;
        sel.bind()?.save(concat!(env!("OUT_DIR"), "/f.pdb"))?;

        // let mut h = FileHandler::create(concat!(env!("OUT_DIR"), "/f.pdb"))?;
        // h.write(&sel)?;
         
         
        println!("test {:p}",&sel);
        println!("out {}",env!("OUT_DIR"));

        //h.write(&sel)?;
        Ok(())
    }

    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let sel = make_sel_prot()?;
        sel.bind_mut()?.unwrap_connectivity_dim(0.2, PBC_FULL)?;

        let mut h = FileHandler::create(concat!(env!("OUT_DIR"), "/unwrapped.pdb"))?;
        h.write_topology(&sel.bind()?)?;
        h.write_state(&sel.bind()?)?;
        Ok(())
    }

    #[test]
    fn eigen_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let sel2 = make_sel_prot()?;
        
        sel2.bind_mut()?.rotate(&Vector3f::x_axis(), 80.0_f32.to_radians());

        sel2.bind()?.save(concat!(env!("OUT_DIR"), "/sel2.pdb"))?;
        sel1.bind()?.save(concat!(env!("OUT_DIR"), "/sel1_before.pdb"))?;
        println!("Initial RMSD:{}", rmsd_mw(&sel1.bind()?, &sel2.bind()?)?);

        let m = fit_transform(&sel1.bind()?, &sel2.bind()?)?;
        println!("{m}");

        sel1.bind_mut()?.apply_transform(&m);

        sel1.bind()?.save(concat!(env!("OUT_DIR"), "/sel1_after.pdb"))?;
        println!("Final RMSD:{}", rmsd_mw(&sel1.bind()?, &sel2.bind()?)?);
        Ok(())
    }

    #[test]
    fn sasa_test() -> anyhow::Result<()> {
        let sel1 = make_sel_all()?;
        let res = sel1.bind()?.sasa();
        println!(
            "Sasa: {a}, Volume: {v}",
            a = res.total_area(),
            v = res.total_volume()
        );
        Ok(())
    }

    #[test]
    fn tets_gyration() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let g = sel1.bind()?.gyration_pbc()?;
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
        let (moments, axes) = sel1.bind()?.inertia_pbc()?;
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
        let tr = sel1.bind()?.principal_transform()?;
        println!("Transform: {tr}");

        let (_, axes) = sel1.bind()?.inertia()?;
        println!("Axes before: {axes}");

        sel1.bind_mut()?.apply_transform(&tr);

        let (_, axes) = sel1.bind()?.inertia()?;
        println!("Axes after: {axes}");

        sel1.bind()?.save(concat!(env!("OUT_DIR"), "/oriented.pdb"))?;

        Ok(())
    }

    #[test]
    fn test_builder_append_from_self() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let n = top.len();
        let builder = System::new(top, st)?;

        let sel = builder.select("resid 550:560")?;
        let added = sel.len();
        builder.append(&sel);
        let all = builder.select_all()?;
        assert_eq!(all.len(), n + added);
        Ok(())
    }

    #[test]
    fn test_select_from_vec() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let src = System::new(top, st)?;
        let last = src.len() - 1;
        let _sel = src.select([1usize, 2, last].as_slice())?;
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_builder_remove_from_self() {
        let (top, st) = FileHandler::open("tests/protein.pdb")
            .unwrap()
            .read()
            .unwrap();
        let n = top.len();
        let builder = System::new(top, st).unwrap();
        let sel = builder.select("resid 550:560").unwrap();
        let removed = sel.len();
        builder.remove(&sel).unwrap(); // Should break here
        let all = builder.select_all().unwrap();
        assert_eq!(all.len(), n - removed);
    }

    #[test]
    #[should_panic]
    fn test_builder_fail_invalid_selection() {
        let (top, st) = FileHandler::open("tests/protein.pdb")
            .unwrap()
            .read()
            .unwrap();
        let builder = System::new(top, st).unwrap();
        let sel = builder.select("resid 809").unwrap(); // last residue
        builder.remove(&sel).unwrap();
        // Trying to call method on invalid selection
        let _cm = sel.center_of_mass();
    }

    // #[test]
    // fn as_parrallel_test() -> anyhow::Result<()> {
    //     let top = TOPST.0.clone();
    //     let st = TOPST.1.clone();
    //     let ser = System::new(top, st)?;
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
        let ser = System::new(top, st)?;
        let ser_sel = ser.select_all()?;

        let mut parts = ser_sel.split_par(|p| Some(p.atom.resindex))?;
        let coms = parts
            .par_iter()
            .map(|sel| sel.center_of_mass())
            .collect::<Result<Vec<_>, _>>()?;

        parts
            .par_iter()
            .enumerate()
            .for_each(|(i, sel)| sel.translate(&(-coms[i].coords)));

        println!("{:?}", coms);

        Ok(())
    }
}
