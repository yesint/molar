use molar::prelude::*;

const PDB: &str = "inp.pdb";
const XTC: &str = "traj_comp.xtc";

fn test1() {
    let t = std::time::Instant::now();

    let src = Source::serial_from_file(PDB).unwrap();
    let ref_sel = src.select_all().unwrap();
    let mut cur_sel = src.select_all().unwrap();

    let mut rmsd = vec![];

    let trj = FileHandler::open(XTC).unwrap().into_iter();
    for st in trj.take(500) {
        cur_sel.set_state(st).unwrap();
        let tr = MeasureMasses::fit_transform(&cur_sel, &ref_sel).unwrap();
        cur_sel.apply_transform(&tr);
        rmsd.push( MeasurePos::rmsd(&cur_sel, &ref_sel).unwrap() );
    }

    println!("Elapsed: {}",t.elapsed().as_secs_f32());
}

fn test2() {
    let t = std::time::Instant::now();

    let src = Source::serial_from_file(PDB).unwrap();
    let mut sel = src.select_str("within 1.0 of protein").unwrap();
    let mut cm = vec![];
    let trj = FileHandler::open(XTC).unwrap().into_iter();
    for st in trj.take(500) {
        sel.set_state(st).unwrap();
        cm.push( sel.center_of_mass().unwrap() );
    }

    println!("Elapsed: {}",t.elapsed().as_secs_f32());
}

fn test3() {
    let t = std::time::Instant::now();

    let src = Source::serial_from_file(PDB).unwrap();
    let mut sel = src.select_str("protein").unwrap();

    let in_trj = FileHandler::open(XTC).unwrap().into_iter();
    let mut out_trj = FileHandler::create("target/.extracted.dcd").unwrap();
    for st in in_trj.take(500) {
        sel.set_state(st).unwrap();
        out_trj.write_state(&sel).unwrap();
    }

    println!("Elapsed: {}",t.elapsed().as_secs_f32());
}

#[cfg(test)]
mod tests {
    use crate::{test1, test2, test3};

    #[ignore]
    #[test]
    fn bench1() {
        test1();
    }

    #[ignore]
    #[test]
    fn bench2() {
        test2();
    }

    #[ignore]
    #[test]
    fn bench3() {
        test3();
    }
}