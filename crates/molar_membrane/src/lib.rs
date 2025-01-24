use std::{ops::Range, sync::Arc};

use molar::prelude::*;

struct LipidSpecies {
    name: String,
    whole_sel_str: String,
    head_marker_subsel_str: String,
    mid_marker_subsel_str: String,
    tails: Vec<TailDescr>,
    head_marker_range: Range<usize>,
    mid_marker_range: Range<usize>,
}

struct TailDescr {
    def_str: String,
    offsets: Vec<usize>,
}

struct LipidMolecule {    
    species: Arc<LipidSpecies>,    
    head_marker: Pos,
    mid_marker: Pos,
    tail_marker: Pos,
    tails: Vec<Pos>,
}

impl LipidMolecule {
    fn update<K: UserCreatableKind>(&mut self, sel: Sel<K>) {
        sel.subsel_local_range(self.species.head_marker_range.clone());
    }
}


