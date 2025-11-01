use crate::prelude::*;

pub trait SelectionLike {
    fn len(&self) -> usize;
    fn get_unchecked(&self,i: usize) -> usize;
    fn iter(&self) -> impl Iterator<Item=usize>;
}

pub struct System {
    top: Topology,
    st: State,
}

impl System {
    pub fn len(&self) -> usize {
        self.st.coords.len()
    }

    fn check_sel(&self, sel: &impl SelectionLike) -> Result<(), SelectionError> {
        if sel.len() == 0 {
            Err(SelectionError::EmptySlice)
        } else if sel.get_unchecked(sel.len()-1) >= self.len() {
            Err(SelectionError::EmptySlice)
        } else {
            Ok(())
        }
    }

    pub fn iter_pos<'a>(&'a self, sel: &'a impl SelectionLike) ->Result<impl Iterator<Item=&'a Pos>,SelectionError> {
        self.check_sel(sel)?;
        Ok(sel.iter().map(|i| unsafe{self.st.coords.get_unchecked(i)}))
    }

    //pub fn get_pos(&self, i: usize)   
}

#[cfg(test)]
mod tests {
    #[test]
    fn test1() {
        let sys = System::default();
        let sel = sys.select("name CA");
        let cm = sys.com(&sel);
        sys.translate(&sel, vec);
        // 1
        let at1 = sys.get_atom(sel.first());
        let at2 = sys.get_atom(sel.get(10));
        let at3 = sys.get_atom(sel.get(20));
        // 2
        sys.with(&sel, |bound|{
            at1 = bound.get(0);
            at2 = bound.get(1);
            at3 = bound.get(2);
            let cm = bound.com();
        });
        // 3
        let b = sys.bind(&sel);
        at1 = b.get(0);
        at2 = b.get(1);
        at3 = b.get(2);
        
    }
}