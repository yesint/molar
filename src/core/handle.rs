use std::sync::{Arc,RwLock,RwLockReadGuard,RwLockWriteGuard};
use std::{str::FromStr, borrow::{Borrow, BorrowMut}, cell::{RefCell,Ref,RefMut}, ops::Deref, ops::DerefMut, rc::Rc};

/// A trait for universal reference-counted lock with internal mutability 
/// allowing multiple readers or a single writer, 
/// which may represent either `Rc<RefCell<T>>` or `Arc<RwLock<T>>`. 
/// 
/// # Basics
/// The `Rc<RefCell<T>>` and `Arc<RwLock<T>>` are idiomatic ways of expressing
/// internal mutability in Rust in the single-threaded and multi-threaded
/// scenarious respectively. Both types provide essentially the same interface
/// with multiple readers or a single writer for their inner data.
/// 
/// However, despite their similarity, these types do not share any common trait
/// and thus can't be used in a generic way: one can't define a function argument 
/// or a struct field that could be either `Rc<RefCell<T>>` or `Arc<RwLock<T>>`.
/// 
/// `UniRcLock` provides a common trait 
/// for `Rc<RefCell<T>>` and `Arc<RwLock<T>>` so that they could be used in
/// a generic way in single-threaded and multi-threaded scenarious alike.
/// 
/// # Performance
/// `UniRcLock` is a zero-cost abstraction.
/// 
/// # Limitations
/// An ability to recover from lock poisoning in `RwLock<T>` is lost
/// when using `UniRcLock`. The methods `read()` and `write()` will panic if 
/// the lock is poisoned.

pub trait UniRcLock<'a,T> {
    type OutRead: Deref<Target=T>+'a;
    type OutWrite: DerefMut<Target=T>+'a;
    /// Obtain a scoped guard for reading
    fn read(&'a self) -> Self::OutRead;
    /// Obtain a scoped guard for writing
    fn write(&'a self) -> Self::OutWrite;
}

// Implementation for Rc<RefCell<T>>
impl<'a,T: 'a> UniRcLock<'a,T> for Rc<RefCell<T>> {
    type OutRead = Ref<'a,T>;
    type OutWrite = RefMut<'a,T>;
    
    fn read(&'a self) -> Self::OutRead {
        Rc::deref(self).borrow()
    }
    
    fn write(&'a self) -> Self::OutWrite {
        Rc::deref(self).borrow_mut()
    }
}

// Implementation for Arc<RwLock<T>>
impl<'a,T: 'a> UniRcLock<'a,T> for Arc<RwLock<T>> {
    type OutRead = RwLockReadGuard<'a,T>;
    type OutWrite = RwLockWriteGuard<'a,T>;
    
    fn read(&'a self) -> Self::OutRead {
        Arc::deref(self).read().expect("Read lock should not be poisoned")
    }
    
    fn write(&'a self) -> Self::OutWrite {
        Arc::deref(self).write().expect("Write lock should not be poisoned")
    }
}

#[cfg(test)]
mod tests {
    use std::{rc::Rc, cell::{RefCell, Ref}, sync::{Arc, RwLock}};

    use super::UniRcLock;

    #[derive(Debug)]
    struct State {val: i32}
    //-----------------
    // Generic handler
    //-----------------
    #[derive(Debug,Clone)]
    struct StateHandler<T: for<'a> UniRcLock<'a,State>> {
        state: T,
    }
    
    impl<T> StateHandler<T> 
    where T: for<'a> UniRcLock<'a, State>
    {
        fn new(val: T) -> Self {
            Self{state: val}
        }        
    }

    #[test]
    fn rc() {
        let st1 = Rc::new(RefCell::new(State{val: 42}));
        st1.write().val +=1;
        println!("{:?}", st1.read());
    }
    #[test]
    fn arc() {
        let st2 = Arc::new(RwLock::new(State{val: 42}));
        st2.write().val +=1;
        println!("{:?}", st2.read());
    }

    #[test]
    fn rc_handler() {
        let st1 = Rc::new(RefCell::new(State{val: 42}));
        let sth1 = StateHandler::new(st1);
        println!("{:?}", sth1.state);
    }

    #[test]
    fn arc_handler() {
        let st2 = Arc::new(RwLock::new(State{val: 42}));
        let sth2 = StateHandler::new(st2);
        sth2.state.write().val += 1;
        println!("{:?}", sth2.state);
    }

    #[test]
    fn threads_test_arc() {
        use std::thread;
        let st2 = Arc::new(RwLock::new(State{val: 42}));
        let sth2 = StateHandler::new(st2);
        
        let threads: Vec<_> = (0..10).map(|i| {
            let h = sth2.clone();
            thread::spawn(move || {
                h.state.write().val += 1;
                println!("Thread #{i} incremented");
            })
        }).collect();

        for t in threads {
            t.join().unwrap();
        }

        println!("Result: {}", sth2.state.read().val);
    }

}