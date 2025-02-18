fn main() {
    // We need to pass different VMDPLUGIN flag for each plugin 
    // cc don't allow per-file defines, so we instead compile each
    // plugin into an object file and them add these objects to a final lib
    let ofiles: Vec<_> = ["pdb","dcd","xyz"].into_iter().map(|pl|
        cc::Build::new()
        .file(format!("molfile/{pl}plugin.c"))
        .define("VMDPLUGIN",format!("{pl}plugin").as_str())
        .include("molfile")
        .pic(true)
        .warnings(false)
        .compile_intermediates()
    ).flatten().collect();
    
    // No compile the whole library
    cc::Build::new()
    .define("VMDUSECONECTRECORDS", None)
    .file("molfile/register_plugins.c")
    .objects(ofiles) // Add precompiled object files
    .include("molfile")
    .pic(true)
    .warnings(false)
    .compile("molfile");
}