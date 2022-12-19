use pteros::MolFileHandler;

fn main() {
    let mut h = MolFileHandler::new("colored.pdb");
    let f = h.read_structure().unwrap();
    //for a in f.atoms {
    //    println!("{:?}",a);
    //}
    let st = h.read_state().unwrap(); 
    for c in st.coords {
        println!("{:?}",c);
    }
}
