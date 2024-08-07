mod kinds;
mod utils;
mod source;
mod source_parallel;
mod sel;
mod sel_split;

pub use kinds::*;
pub use source::*;
pub use source_parallel::*;
pub use sel::*;
pub use sel_split::*;
pub(crate) use utils::*;

use thiserror::Error;

use crate::io::FileIoError;

use super::{selection_parser::SelectionParserError, BuilderError};

#[derive(Error,Debug)]
#[error("topology and state have different sizes ({0},{1})")]
pub struct TopologyStateSizes(usize,usize);

#[derive(Error,Debug)]
pub enum SelectionError {
    #[error("selection parser failed")]
    Parser(#[from] SelectionParserError),
    
    #[error(transparent)]
    DifferentSizes(#[from] TopologyStateSizes),
    
    #[error("creating selection from expression {expr_str}")]
    FromExpr{
        expr_str: String,
        source: SelectionIndexError,
    },
    
    #[error("creating selection from range {first}:{last}")]
    FromRange{
        first: usize,
        last: usize,
        source: SelectionIndexError,
    },

    #[error("invalid local sub-range: {0}:{1}, valid range: 0:{2}")]
    LocalRange(usize,usize,usize),
    
    #[error("creating selection from vec {first}..{last} of size {size}")]
    FromVec{
        first: usize,
        last: usize,
        size: usize,
        source: SelectionIndexError
    },

    #[error("pbc operation on selection withon periodic box")]
    NoPbc,

    #[error("index {0} is beyond the allowed range 0:{1}")]
    OutOfBounds(usize, usize),

    #[error("selection index {0}:{1} is outside the source range: 0:{2}")]
    IndexCheck(usize, usize, usize),

    #[error("index {0} is already used by other selection")]
    OverlapCheck(usize),

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
}

#[derive(Error,Debug)]
pub enum SelectionIndexError {
    #[error("selection index is empty")]
    IndexEmpty,
    #[error("selection indeces {0}:{1} are out of allowed range 0:{2}")]
    IndexOutOfBounds(usize,usize,usize),
}

//############################################################
//#  Tests
//############################################################

#[cfg(test)]
mod tests {    
    use crate::prelude::*;    
    pub use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
    use super::*;

    pub fn read_test_pdb() -> (Topology, State) {
        FileHandler::open("tests/protein.pdb").unwrap().read().unwrap()
    }

    #[test]
    fn builder_overlap() -> anyhow::Result<()> {
        let (top, st) = read_test_pdb();
        let b = Source::new(top, st)?;
        // Create two overlapping selections
        let _sel1 = b.select_from_iter(0..10)?;
        let _sel2 = b.select_from_iter(5..15)?;
        Ok(())
    }

    #[test]
    fn builder_par_no_overlap() {
        let (top, st) = read_test_pdb();
        let b = Source::new(top, st).unwrap();
        // Create two non-overlapping selections.
        let _sel1 = b.select_from_iter(0..10).unwrap();
        let _sel2 = b.select_from_iter(11..15).unwrap();
    }

    #[test]
    #[should_panic]
    fn builder_par_overlap() {
        let (top, st) = read_test_pdb();
        let mut b = SourceParallel::new_mut(top, st).unwrap();
        // Create two overlapping selections. This must fail!
        let _sel1 = b.add_from_iter(0..10).unwrap();
        let _sel2 = b.add_from_iter(5..15).unwrap();
    }

    #[test]
    fn builder_par_test() -> anyhow::Result<()> {
        let (top, st) = read_test_pdb();
        let mut b = SourceParallel::new_mut(top, st)?;
        // Create two valid non-overlapping selections.
        b.add_from_iter(0..10).unwrap();
        b.add_from_iter(11..15).unwrap();
        b.add_from_iter(15..25).unwrap();
        // Process them
        let v = Vector3f::new(1.0, 2.0, 3.0);
        let res = b.par_iter().map(|sel| {
            sel.translate(&v);
            Ok(sel.center_of_mass()?)
        }).collect::<anyhow::Result<Vec<_>>>()?;
        
        println!("cm1 = {:?}",res);
        Ok(())
    }

    #[test]
    fn builder_par_to_serial() -> anyhow::Result<()> {
        let (top, st) = read_test_pdb();
        let mut b = SourceParallel::new_mut(top, st)?;
        // Create two valid non-overlapping selections.
        for i in 0..30 {
            b.add_from_iter(10*i..10*(i+1)).unwrap();
        }
        // Process them
        let v = Vector3f::new(1.0, 2.0, 3.0);
        let mut res: Vec<_> = b.collect_par(|sel| {
            let cm = sel.center_of_mass()?;
            sel.translate(&cm.coords);
            println!("thread: {}", rayon::current_thread_index().unwrap());
            Ok::<_,MeasureError>(sel.center_of_mass()?)
        })?;
        
        println!("cm before = {:?}",res);
        
        let (_,sels) = b.into_serial_with_sels();
        sels[0].translate(&v);
        res[0] = sels[0].center_of_mass()?;
        
        println!("cm after = {:?}",res);

        Ok(())
    }

    fn make_sel_all() -> anyhow::Result<Sel<MutableSerial>> {
        let (top, st) = read_test_pdb();
        let b = Source::new(top, st)?;
        let sel = b.select_all()?;
        Ok(sel)
    }

    fn make_sel_prot() -> anyhow::Result<Sel<MutableSerial>> {
        let (top, st) = read_test_pdb();
        let b = Source::new(top, st)?;
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
    fn test_searcher_single() {
        let sel = make_sel_prot().unwrap();
        let mut searcher = DistanceSearcherSingle::new_periodic(
            0.3, 
            sel.iter_particle().map(|p| (p.id,*p.pos)), 
            sel.get_box().unwrap(), 
            &PBC_FULL
        );

        // Enforce parallel
        searcher.set_serial_limit(0);
        let v1: Vec<usize> = searcher.search();

        // Enforce serial
        searcher.set_serial_limit(1e10 as usize);
        let v2: Vec<usize> = searcher.search();

        assert_eq!(v1.len(),v2.len());
    }

    #[test]
    fn test_searcher_single_vdw() {
        let sel = make_sel_prot().unwrap();
        let searcher = DistanceSearcherSingle::new_vdw_periodic(
            sel.iter_particle().map(|p| (p.id,*p.pos)), 
            sel.iter_particle().map(|p| p.atom.vdw()),
            sel.get_box().unwrap(), 
            &PBC_FULL
        );

        let v1: Vec<usize> = searcher.search();

        assert!(v1.len()>0);
    }


    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let sel = make_sel_prot()?;
        sel.unwrap_connectivity_dim(0.2, &PBC_FULL)?;

        let mut h = FileHandler::create(concat!(env!("OUT_DIR"),"/unwrapped.pdb"))?;
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

        sel1.save(concat!(env!("OUT_DIR"),"/sel1_after.pdb"))?;
        println!("Final RMSD:{}", Sel::rmsd_mw(&sel1, &sel2)?);
        Ok(())
    }

    #[test]
    fn split_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        for res in sel1.split_contig(|p| p.atom.resid) {
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

        sel1.save(concat!(env!("OUT_DIR"),"/oriented.pdb"))?;

        Ok(())
    }    

    #[test]
    #[should_panic]
    fn fail_on_empty_selection() {
        let (top, st) = FileHandler::open("tests/protein.pdb").unwrap().read().unwrap();
        let mut src = SourceParallel::new_mut(top, st).unwrap();
        src.add_str("resid 5").unwrap();        
    }

    #[test]
    fn test_builder_append_from_self() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let n = top.num_atoms();
        let mut builder = Source::new_builder(top, st)?;
        let sel = builder.select_str("resid 550:560")?;
        let added = sel.len();
        builder.append(&sel);
        let all = builder.select_all()?;
        assert_eq!(all.len(),n+added);
        Ok(())
    }

    #[test]
    fn test_builder_remove_from_self() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let n = top.num_atoms();
        let mut builder = Source::new_builder(top, st)?;
        let sel = builder.select_str("resid 550:560")?;
        let removed = sel.len();
        builder.remove(&sel)?;
        let all = builder.select_all()?;
        assert_eq!(all.len(),n-removed);
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_builder_fail_invalid_selection() {
        let (top, st) = FileHandler::open("tests/protein.pdb").unwrap().read().unwrap();
        let mut builder = Source::new_builder(top, st).unwrap();
        let sel = builder.select_str("resid 809").unwrap(); // last residue
        builder.remove(&sel).unwrap();
        // Trying to call method on invalid selection
        let _cm = sel.center_of_mass();
    }

}
