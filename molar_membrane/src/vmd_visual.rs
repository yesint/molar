use molar::prelude::*;
use std::path::Path;

pub struct VmdVisual {
    buf: String,
}

impl VmdVisual {
    pub fn new() -> Self {
        Self { buf: String::new() }
    }

    pub fn save_to_file(&self, fname: impl AsRef<Path>) -> anyhow::Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(fname)?;
        writeln!(f, "{}", self.buf)?;
        Ok(())
    }

    pub fn sphere(&mut self, point: &Pos, radius: f32, color: &str) {
        use std::fmt::Write;
        writeln!(self.buf, "draw color {color}").unwrap();
        writeln!(
            self.buf,
            "draw sphere \"{} {} {}\" radius {radius} resolution 12",
            point.x * 10.0,
            point.y * 10.0,
            point.z * 10.0
        )
        .unwrap();
    }

    pub fn arrow(&mut self, point: &Pos, dir: &Vector3f, color: &str) {
        use std::fmt::Write;

        const LENGTH: f32 = 5.0;

        let p1 = point * 10.0;
        let p2 = p1 + dir * 0.5 * LENGTH;
        let p3 = p1 + dir * 0.7 * LENGTH;

        writeln!(self.buf, "draw color {color}").unwrap();

        writeln!(
            self.buf,
            "draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2 resolution 12",
            p1.x, p1.y, p1.z, p2.x, p2.y, p2.z
        )
        .unwrap();

        writeln!(
            self.buf,
            "draw cone \"{} {} {}\" \"{} {} {}\" radius 0.4 resolution 12\n",
            p2.x, p2.y, p2.z, p3.x, p3.y, p3.z
        )
        .unwrap();
    }

    pub fn cylinder(&mut self, point1: &Pos, point2: &Pos, color: &str) {
        use std::fmt::Write;

        let p1 = point1 * 10.0;
        let p2 = point2 * 10.0;

        writeln!(self.buf, "draw color {color}").unwrap();

        writeln!(
            self.buf,
            "draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2 resolution 12",
            p1.x, p1.y, p1.z, p2.x, p2.y, p2.z
        )
        .unwrap();
    }
}
