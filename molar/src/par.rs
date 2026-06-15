//! Parallelism shim.
//!
//! Native targets re-export rayon's iterator traits unchanged. `wasm32` has no
//! threads without `SharedArrayBuffer`, so molar runs **single-threaded** there:
//! the "parallel" entry points fall back to plain std iterators, mirroring only
//! the rayon surface molar actually uses (`par_iter` / `into_par_iter` /
//! `with_min_len` / `from_par_iter`). All downstream combinators (`map`, `collect`,
//! `cloned`, `enumerate`, `zip`, …) are ordinary `Iterator` methods, so the call
//! sites compile unchanged under both targets.

#[cfg(not(target_arch = "wasm32"))]
pub use rayon::iter::{
    FromParallelIterator, IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    ParallelIterator,
};

#[cfg(target_arch = "wasm32")]
pub use serial::*;

#[cfg(target_arch = "wasm32")]
mod serial {
    // A "parallel" iterator is just a std iterator on wasm, so map/collect/cloned/
    // enumerate/zip/filter resolve to std `Iterator`. Aliasing the trait names to it
    // both brings those methods into scope at the call sites and keeps
    // `impl IndexedParallelIterator<Item = T>` return types valid.
    pub use std::iter::Iterator as IndexedParallelIterator;
    pub use std::iter::Iterator as ParallelIterator;

    /// `slice.par_iter()` / `vec.par_iter()` → `.iter()`.
    pub trait IntoParallelRefIterator<'data> {
        type Item;
        type Iter: Iterator<Item = Self::Item>;
        fn par_iter(&'data self) -> Self::Iter;
    }
    impl<'data, T: 'data> IntoParallelRefIterator<'data> for [T] {
        type Item = &'data T;
        type Iter = std::slice::Iter<'data, T>;
        fn par_iter(&'data self) -> Self::Iter {
            self.iter()
        }
    }
    impl<'data, T: 'data> IntoParallelRefIterator<'data> for Vec<T> {
        type Item = &'data T;
        type Iter = std::slice::Iter<'data, T>;
        fn par_iter(&'data self) -> Self::Iter {
            self.iter()
        }
    }

    /// `x.into_par_iter()` → `x.into_iter()`.
    pub trait IntoParallelIterator {
        type Item;
        type Iter: Iterator<Item = Self::Item>;
        fn into_par_iter(self) -> Self::Iter;
    }
    impl<T: IntoIterator> IntoParallelIterator for T {
        type Item = T::Item;
        type Iter = T::IntoIter;
        fn into_par_iter(self) -> Self::Iter {
            self.into_iter()
        }
    }

    /// `iter.with_min_len(n)` → no-op (a serial iterator has no work-splitting).
    pub trait WithMinLen: Iterator + Sized {
        fn with_min_len(self, _min: usize) -> Self {
            self
        }
    }
    impl<I: Iterator> WithMinLen for I {}

    /// `C::from_par_iter(it)` → `it.into_iter().collect()`; a blanket over every
    /// `FromIterator` target so `C: FromParallelIterator<T>` bounds stay satisfiable.
    pub trait FromParallelIterator<T>: Sized {
        fn from_par_iter<I: IntoIterator<Item = T>>(it: I) -> Self;
    }
    impl<T, C: FromIterator<T>> FromParallelIterator<T> for C {
        fn from_par_iter<I: IntoIterator<Item = T>>(it: I) -> Self {
            it.into_iter().collect()
        }
    }
}
