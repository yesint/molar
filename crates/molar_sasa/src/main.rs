#[cxx::bridge]
mod ffi {
    unsafe extern "C++" {
        include!("molar_sasa/powersasa/wrapper.h");

        type PowerSasaWrapper;

        fn new_power_sasa(
            crd_prov: &mut RustCoordProvider,
            vdw_prov: &mut RustVdwProvider,
            sz: usize,
        ) -> UniquePtr<PowerSasaWrapper>;

        fn get_next_coord(&self);
    }

    extern "Rust" {
        type RustCoordProvider;
        fn get(self: &mut RustCoordProvider) -> *mut f32;
        fn advance(self: &mut RustCoordProvider);

        type RustVdwProvider;
        fn get(self: &RustVdwProvider) -> f32;
        fn advance(self: &mut RustVdwProvider);
    }
}

struct RustCoordProvider {
    val: [f32; 3],
}

impl RustCoordProvider {
    pub fn get(&mut self) -> *mut f32 {
        self.val.as_mut_ptr()
    }

    pub fn advance(&mut self) {
        self.val[0] += 1.0;
    }
}
struct RustVdwProvider {
    val: f32,
}

impl RustVdwProvider {
    pub fn get(&self) -> f32 {
        self.val
    }

    pub fn advance(&mut self) {
        self.val += 1.0;
    }
}

fn main() {
    let mut crd_prov = RustCoordProvider {
        val: [0.0, 0.0, 0.0],
    };
    let mut vdw_prov = RustVdwProvider { val: 42.0 };
    let h = ffi::new_power_sasa(&mut crd_prov, &mut vdw_prov, 3);
    h.get_next_coord();
    h.get_next_coord();
    h.get_next_coord();
}
