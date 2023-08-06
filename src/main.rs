use pteros::io::*;
//use pteros::{MolfileSingleFrame, MolfileStructure};

fn test_pdb() {
    let mut h = VmdMolFileHandler::new("colored.pdb");
    let f = h.read_structure().unwrap();
    for _a in f.atoms {
        //println!("{:?}",a);
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
    let mut h = XtcFileHandler::new("no_ATP.xtc");
    h.open_read().unwrap();
    let st = h.read_next_state().unwrap().unwrap();
    println!("{}", st.coords.len());
    println!("{:?}",st.box_);
    
    let mut wh = XtcFileHandler::new("new.xtc");
    wh.open_write().unwrap();
    for i in 0..10 {
        wh.write_next_state(&st);
    }
}

fn main() {
    test_xtc();
}
