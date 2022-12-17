use pteros::PdbFile;

fn main() {
    let mut pdb = PdbFile::new().open_read("colored.pdb").read();
}
