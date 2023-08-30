use pteros::io::*;
//use pteros::{MolfileSingleFrame, MolfileStructure};

#[test]
fn test_pdb() {
    let mut h = VmdMolFileHandler::new_reader("colored.pdb").unwrap();
    let st = h.read_structure().unwrap();
    for _a in st.atoms {
        println!("{:?}",_a);
    }

    for _i in 0..2 {
        let st = h.read_state().unwrap();
        
        //if let Some(s) = st {
            println!("{}", st.coords.len());
            for c in st.coords {
                println!("{:?}",c);
            }
        //} else {
        //    println!("EOF reached");
        //}
    }
}

#[test]
fn test_xtc() {
    let mut h = XtcFileHandler::new_reader("no_ATP.xtc").unwrap();
    let st = h.read_next_state().unwrap().unwrap();
    println!("{}", st.coords.len());
    println!("{:?}",st.box_);
    
    let mut wh = XtcFileHandler::new_writer("new.xtc").unwrap();
    for _ in 0..10 {
        wh.write_next_state(&st).unwrap();
    }
}

fn main() {
    let h = FileHandler::new_reader("colored.pdb").unwrap();
    if let FileHandler::Pdb(mut e) = h {
        let st = e.read_structure().unwrap();
        println!("{}",st.atoms.len());
    }

}
