use crate::prelude::*;
use rayon::iter::{FromParallelIterator, IntoParallelIterator};

#[derive(Debug, Default)]
pub struct SearchConnectivity(rustc_hash::FxHashMap<usize, Vec<usize>>);

impl SearchConnectivity {
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> SearchConnectivityIter {
        SearchConnectivityIter(self.0.iter())
    }

    fn from_iter(iter: impl IntoIterator<Item = (usize, usize)>) -> Self {
        let mut res = Self::default();
        for (i, j) in iter {
            if res.0.contains_key(&i) {
                res.0.get_mut(&i).unwrap().push(j);
            } else {
                res.0.insert(i, vec![j]);
            }

            if res.0.contains_key(&j) {
                res.0.get_mut(&j).unwrap().push(i);
            } else {
                res.0.insert(j, vec![i]);
            }
        }
        res
    }
}

impl FromIterator<(usize, usize)> for SearchConnectivity {
    fn from_iter<T: IntoIterator<Item = (usize, usize)>>(iter: T) -> Self {
        Self::from_iter(iter)
    }
}

impl FromParallelIterator<(usize, usize)> for SearchConnectivity {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = (usize, usize)>,
    {
        let v = par_iter.into_par_iter().collect::<Vec<_>>();
        Self::from_iter(v.iter().cloned())
    }
}

pub struct SearchConnectivityIter<'a>(std::collections::hash_map::Iter<'a, usize, Vec<usize>>);

impl<'a> Iterator for SearchConnectivityIter<'a> {
    type Item = (&'a usize, &'a Vec<usize>);
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}

impl IntoIterator for SearchConnectivity {
    type Item = (usize, Vec<usize>);
    type IntoIter = std::collections::hash_map::IntoIter<usize, Vec<usize>>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl std::ops::Index<usize> for SearchConnectivity {
    type Output = Vec<usize>;
    fn index(&self, i: usize) -> &Self::Output {
        &self.0[&i]
    }
}

//---------------------------------------------------------------
#[derive(Debug, Default)]
pub struct LocalConnectivity(pub Vec<Vec<usize>>);

impl LocalConnectivity {
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn from_iter(iter: impl IntoIterator<Item = (usize, usize)>, table_size: usize) -> Self {
        let mut res = Self(vec![vec!(); table_size]);
        // Now populate
        for (i, j) in iter {
            res.0[i].push(j);
            res.0[j].push(i);
        }
        res
    }
}
