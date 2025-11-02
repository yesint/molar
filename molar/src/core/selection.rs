mod selection_def;
mod utils;

use std::num::ParseIntError;

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
    #[error(transparent)]
    Bind(#[from] BindError),

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

    fn make_sel_all() -> anyhow::Result<(System,Sel)> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = System::new(top, st)?;
        let sel = b.select_all();
        Ok((b,sel))
    }

    fn make_sel_prot() -> anyhow::Result<(System,Sel)> {
        let top = TOPST.0.clone();
        let st = TOPST.1.clone();
        let b = System::new(top, st)?;
        let sel = b.select("not resname TIP3 POT CLA")?;
        Ok((b,sel))
    }

    #[test]
    fn test_measure() -> anyhow::Result<()> {
        let (sys,sel) = make_sel_all()?;
        println!("before {}", sel.bind(&sys)?.iter_pos().next().unwrap());

        let (minv, maxv) = sel.bind(&sys)?.min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}", sel.bind(&sys)?.iter_pos().next().unwrap());
        println!("{:?}", sel.bind(&sys)?.min_max());
        Ok(())
    }

    #[test]
    fn test_measure_pbc() -> anyhow::Result<()> {
        let (sys,sel) = make_sel_all()?;

        let cm = sel.bind(&sys)?.center_of_mass()?;
        println!("{cm}");
        Ok(())
    }

    #[test]
    fn test_translate() -> anyhow::Result<()> {
        let (mut sys,sel) = make_sel_all()?;

        println!("before {}", sel.bind(&sys)?.iter_pos().next().unwrap());
        sel.bind_mut(&mut sys)?.translate(&Vector3f::new(10.0, 10.0, 10.0));
        println!("after {}", sel.bind(&sys)?.iter_pos().next().unwrap());
        Ok(())
    }

    #[test]
    fn test_write_to_file() -> anyhow::Result<()> {
        let sys = System::from_file("tests/protein.pdb")?;
        let sel = sys.select("name CA")?;
        sel.bind(&sys)?.save(concat!(env!("OUT_DIR"), "/f.pdb"))?;

        // let mut h = FileHandler::create(concat!(env!("OUT_DIR"), "/f.pdb"))?;
        // h.write(&sel)?;
         
         
        println!("test {:p}",&sel);
        println!("out {}",env!("OUT_DIR"));

        //h.write(&sel)?;
        Ok(())
    }

    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let (mut sys,sel) = make_sel_prot()?;
        sys.unwrap_connectivity_dim(&sel,0.2, PBC_FULL)?;

        let mut h = FileHandler::create(concat!(env!("OUT_DIR"), "/unwrapped.pdb"))?;
        h.write_topology(&sel.bind(&sys)?)?;
        h.write_state(&sel.bind(&sys)?)?;
        Ok(())
    }

    #[test]
    fn eigen_test() -> anyhow::Result<()> {
        let (mut sys,sel1) = make_sel_prot()?;
        let (_,sel2) = make_sel_prot()?;

        sel2.bind_mut(&mut sys)?.rotate(&Vector3f::x_axis(), 80.0_f32.to_radians());

        sys.save_sel(&sel1,concat!(env!("OUT_DIR"), "/sel2.pdb"))?;
        sys.save_sel(&sel2,concat!(env!("OUT_DIR"), "/sel1_before.pdb"))?;
        println!("Initial RMSD:{}", sys.rmsd_mw(&sel1, &sel2)?);

        let m = sys.fit_transform(&sel1, &sel2)?;
        println!("{m}");

        sys.apply_transform(&sel1,&m)?;

        sys.save_sel(&sel1,concat!(env!("OUT_DIR"), "/sel1_after.pdb"))?;
        println!("Final RMSD:{}", sys.rmsd_mw(&sel1, &sel2)?);
        Ok(())
    }

    #[test]
    fn sasa_test() -> anyhow::Result<()> {
        let (sys,sel1) = make_sel_all()?;
        let res = sys.sasa(&sel1)?;
        println!(
            "Sasa: {a}, Volume: {v}",
            a = res.total_area(),
            v = res.total_volume()
        );
        Ok(())
    }

    #[test]
    fn tets_gyration() -> anyhow::Result<()> {
        let (sys,sel1) = make_sel_prot()?;
        let g = sys.gyration_pbc(&sel1)?;
        println!("Gyration radius: {g}");
        Ok(())
    }

    #[test]
    fn test_inertia() -> anyhow::Result<()> {
        let (sys,sel1) = make_sel_all()?;
        //sel1.rotate(
        //    &UnitVector3::new_normalize(Vector3f::new(1.0,2.0,3.0)),
        //    0.45
        //);
        let (moments, axes) = sys.inertia_pbc(&sel1)?;
        println!("Inertia moments: {moments}");
        println!("Inertia axes: {axes:?}");
        Ok(())
        // {0.7308828830718994 0.3332606256008148 -0.5956068634986877}
        // {0.6804488301277161 -0.28815552592277527 0.6737624406814575}
        // {-0.05291106179356575 0.8977214694023132 0.4373747706413269}
    }

    #[test]
    fn test_principal_transform() -> anyhow::Result<()> {
        let (mut sys,sel1) = make_sel_prot()?;
        let tr = sys.principal_transform(&sel1)?;
        println!("Transform: {tr}");

        let (_, axes) = sys.inertia(&sel1)?;
        println!("Axes before: {axes}");

        sys.apply_transform(&sel1,&tr)?;

        let (_, axes) = sys.inertia(&sel1)?;
        println!("Axes after: {axes}");

        sys.save_sel(&sel1,concat!(env!("OUT_DIR"), "/oriented.pdb"))?;

        Ok(())
    }

    #[test]
    fn test_builder_append_from_self() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let n = top.len();
        let mut builder = System::new(top, st)?;

        let sel = builder.select("resid 550:560")?;
        let added = sel.len();
        builder.append_sel(&sel);
        let all = builder.select_all();
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

    // #[test]
    // #[should_panic]
    // fn test_builder_remove_from_self() {
    //     let (top, st) = FileHandler::open("tests/protein.pdb")
    //         .unwrap()
    //         .read()
    //         .unwrap();
    //     let n = top.len();
    //     let builder = System::new(top, st).unwrap();
    //     let sel = builder.select("resid 550:560").unwrap();
    //     let removed = sel.len();
    //     builder.remove(&sel).unwrap(); // Should break here
    //     let all = builder.select_all();
    //     assert_eq!(all.len(), n - removed);
    // }

    // #[test]
    // #[should_panic]
    // fn test_builder_fail_invalid_selection() {
    //     let (top, st) = FileHandler::open("tests/protein.pdb")
    //         .unwrap()
    //         .read()
    //         .unwrap();
    //     let builder = System::new(top, st).unwrap();
    //     let sel = builder.select("resid 809").unwrap(); // last residue
    //     builder.remove(&sel).unwrap();
    //     // Trying to call method on invalid selection
    //     let _cm = sel.center_of_mass();
    // }

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
        let mut ser = System::new(top, st)?;

        let parts = ser.split_par(|p| Some(p.atom.resindex))?;
        let coms = parts.bind_mut(&mut ser)?
            .iter_par()
            .map(|sel| sel.center_of_mass())
            .collect::<Result<Vec<_>, _>>()?;

        parts.bind_mut(&mut ser)?
            .iter_par_mut()
            .enumerate()
            .for_each(|(i, sel)| sel.translate(&(-coms[i].coords)));

        println!("{:?}", coms);

        Ok(())
    }
}
