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
        type RustCoordProvider<'a>;
        fn get(self: &RustCoordProvider) -> *const f32;
        fn advance(self: &mut RustCoordProvider);
        fn reset(self: &mut RustCoordProvider);

        type RustVdwProvider<'a>;
        fn get(self: &RustVdwProvider) -> f32;
        fn advance(self: &mut RustVdwProvider);
        fn reset(self: &mut RustVdwProvider);
    }
}

struct RustCoordProvider<'a> {
    closure: &'a dyn Fn(usize) -> *const f32,
    cur: usize,
}

impl<'a> RustCoordProvider<'a> {
    pub fn new(f: &'a dyn Fn(usize) -> *const f32) -> Self {
        Self { closure: f, cur: 0 }
    }

    pub fn get(&self) -> *const f32 {
        (self.closure)(self.cur)
    }

    pub fn advance(&mut self) {
        self.cur += 1;
    }

    pub fn reset(&mut self) {
        self.cur = 0;
    }
}
struct RustVdwProvider<'a> {
    closure: &'a dyn Fn(usize) -> f32,
    cur: usize,
}

impl<'a> RustVdwProvider<'a> {
    pub fn new(f: &'a dyn Fn(usize) -> f32) -> Self {
        Self { closure: f, cur: 0 }
    }

    pub fn get(&self) -> f32 {
        (self.closure)(self.cur)
    }

    pub fn advance(&mut self) {
        self.cur += 1;
    }

    pub fn reset(&mut self) {
        self.cur = 0;
    }
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test1() {
        let mut vec = vec![1.0f32, 2.0, 3.0];

        let f = |_i| vec.as_ptr();

        let g = |i: usize| i as f32;

        let mut crd_prov = RustCoordProvider::new(&f);
        let mut vdw_prov = RustVdwProvider::new(&g);
        let h = ffi::new_power_sasa(&mut crd_prov, &mut vdw_prov, 3);
        h.get_next_coord();
        h.get_next_coord();
        h.get_next_coord();
    }
}
