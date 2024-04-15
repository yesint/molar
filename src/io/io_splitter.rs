use anyhow::{anyhow, Result};
use crate::core::{StateUArc, TopologyUArc};

pub trait ReadTopAndState {    
    fn read_top_and_state(&mut self) -> Result<(TopologyUArc, StateUArc)>;
    
    //fn write_top_and_state(&mut self, data: &(impl IoIndexProvider + IoTopologyProvider + IoStateProvider)) -> Result<()> {
    //    unimplemented!()
    //}
}

pub struct IoSplitter<T: ReadTopAndState> {
    pub handler: T,
    // For reading
    top: Option<TopologyUArc>,
    state: Option<StateUArc>,
    // For writing
    //top_ref: Option<&'a Topology>,
    //state_ref: Option<&'a State>,
    //index_ref: Option<&'a dyn IndexIterator>,
    // Flag for actual IO operation performed
    io_done: bool,
}

impl<T: ReadTopAndState> IoSplitter<T> {
    pub fn new(handler: T) -> Self {
        IoSplitter {
            handler,
            top: None,
            state: None,
            //top_ref: None,
            //state_ref: None,
            //index_ref: None,
            io_done: false, 
        }
    }

    pub fn read(&mut self) -> Result<(TopologyUArc,StateUArc)> {
        self.io_done = true;
        self.handler.read_top_and_state()
    }

    pub fn read_topology(&mut self) -> Result<TopologyUArc> {
        self.try_read()?;
        self.top.take().ok_or_else(|| anyhow!("Topology already read!"))
    }

    pub fn read_state(&mut self) -> Result<Option<StateUArc>> {
        self.try_read()?;
        Ok(self.state.take())
    }

    fn try_read(&mut self) -> Result<()> {
        if !self.io_done {
            let (top,st) = self.handler.read_top_and_state()?;
            self.top = Some(top);
            self.state = Some(st);
            self.io_done = true;
        }
        Ok(())
    }
    
    /*
    pub fn write(
        &mut self,
        data: &(impl IoIndexProvider + IoTopologyProvider + super::IoStateProvider),
    ) -> Result<()>
    {
        self.handler.write_top_and_state(data)
    }

    
    pub fn write_topology(&mut self, data: &(impl IoIndexProvider+IoTopologyProvider)) -> Result<()> {
        self.top_ref = Some(&data.get_topology());
        if self.index_ref.is_none() {
            self.index_ref = Some(&mut data.get_index());
        }
        Ok(())
    }

    pub fn write_state(&mut self, data: &(impl IoIndexProvider+IoStateProvider)) -> Result<()> {
        self.state_ref = Some(&data.get_state());
        if self.index_ref.is_none() {
            self.index_ref = Some(&mut data.get_index());
        }
        Ok(())
    }


    
    fn try_write(&mut self) -> Result<()>
    {
        if !self.io_done && self.top_ref.is_some() && self.state_ref.is_some() {
            self.handler.write_top_and_state(self)?;
            self.io_done = true;
            Ok(())
        } else {
            bail!("Both Topology and State have to be provided!")
        }
        
    }
    */

}

/*
impl<T: ReadTopAndState> IoIndexProvider for IoSplitter<'_,T> {
    fn get_index(&self) -> impl IndexIterator {
        if let Some(i) = self.index_ref {
            *i
        } else {
            unreachable!()
        }
        
    }
}
*/