use pteros::io::*;
//use pteros::{MolfileSingleFrame, MolfileStructure};

fn test_pdb() {
    let mut h = VmdMolFileHandler::new_reader("colored.pdb");
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

fn test_xtc() {
    let mut h = XtcFileHandler::new_reader("no_ATP.xtc");
    let st = h.read_next_state().unwrap().unwrap();
    println!("{}", st.coords.len());
    println!("{:?}",st.box_);
    
    let mut wh = XtcFileHandler::new_writer("new.xtc");
    for _ in 0..10 {
        wh.write_next_state(&st).unwrap();
    }
}

fn main() {
    test_xtc();
}
