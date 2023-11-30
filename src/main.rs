use molar::io::*;

#[cfg(feature = "gromacs")]
#[test]
fn test_xtc() {
    let mut h = XtcFileHandler::new_reader("tests/no_ATP.xtc").unwrap();
    let st_ref = h.read_next_state().unwrap().unwrap();
    let st = st_ref.read().unwrap();
    println!("{}", st.coords.len());
    println!("{:?}",st.box_);
    
    let mut wh = XtcFileHandler::new_writer("new.xtc").unwrap();
    for _ in 0..10 {
        wh.write_next_state(&st).unwrap();
    }
}

fn main() {
    //let h = TprFileHandler::new_reader("tests/topol.tpr");
}
