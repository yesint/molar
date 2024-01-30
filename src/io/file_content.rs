pub struct FileContent {
    topology: bool,
    single_state: bool,
    trajectory: bool,
    random_access: bool,
    read_only: bool,
}

impl FileContent {
    pub fn new() -> Self {
        FileContent{
            topology: false,
            single_state: false,
            trajectory: false,
            random_access: false,
            read_only: false,
        }
    }

    pub fn with_topology(mut self) -> Self {
        self.topology = true;
        self
    }

    pub fn with_single_state(mut self) -> Self {
        self.single_state = true;
        self
    }

    pub fn with_trajectory(mut self) -> Self {
        self.trajectory = true;
        self
    }

    pub fn with_random_access(mut self) -> Self {
        self.trajectory = true;
        self.random_access = true;
        self
    }

    pub fn with_read_only(mut self) -> Self {
        self.read_only = true;
        self
    }

    pub fn is_topology(&self) -> bool {
        self.topology
    }

    pub fn is_single_state(&self) -> bool {
        self.single_state
    }

    pub fn is_trajectory(&self) -> bool {
        self.trajectory
    }

    pub fn is_random_access(&self) -> bool {
        self.random_access
    }

    pub fn is_read_only(&self) -> bool {
        self.read_only
    }
}
