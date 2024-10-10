fn main() {
    cc::Build::new()
    .cpp(true)
    .file("powersasa/wrapper.cpp")
    .include("powersasa")
    .pic(true)
    .warnings(false)
    .compile("powersasa");
}