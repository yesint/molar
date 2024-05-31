use crate::prelude::*;
use anyhow::Result;

struct Builder {
    source: Source,
}

impl Builder {
    pub fn new(
        topology: TopologyUArc,
        state: StateUArc,
    ) -> anyhow::Result<Self> {
        Ok(Self{source: Source::new(topology,state)?})
    }

    pub fn add<K: SelectionKind>(&self,sel: &Sel<K>) {
        
    }
}