use pteros::MolFileHandler;

fn main() {
    let f = MolFileHandler::new("colored.pdb").read_structure().unwrap();
    for a in f.atoms {
        println!("{:?}",a);
    }
}
