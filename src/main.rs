use pteros::VmdMolFileHandler;
use pteros::{MolfileSingleFrame, MolfileStructure};

fn main() {
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
