fn main() {
    cc::Build::new()
    .file("xdrfile/xdrfile.c")
    .file("xdrfile/xdrfile_xtc.c")
    .file("xdrfile/xdrfile_trr.c")
    .file("xdrfile/xdr_utils.c")
    .file("xdrfile/xdr_seek.c")
    .include("xdrfile")
    .pic(true)
    .warnings(false)
    .compile("xdrfile");
}