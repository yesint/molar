use std::fs::File;
use std::io::Write;

pub type Vector2 = nalgebra::Vector2<f32>;
const TOL: f32 = 1e-10;

struct Vertex {
    pos: Vector2,
    ccw: usize,
    cw: usize,
}

struct CuttingLine {
    pos: Vector2,
    r: f32,
}

impl CuttingLine {
    fn new(point: &Vector2) -> Self {
        let pos = 0.5*point;
        let r = pos.norm();
        Self {pos, r}
    }
}

pub struct VoronoiCell {
    // Vertices
    vert: Vec<Vertex>,
    init_vert: usize,
}

#[derive(Default,Debug)]
struct CutEdge {
    in_ind: usize,
    in_dist: f32,
    out_ind: usize,
    out_dist: f32,
}

impl VoronoiCell {
    pub fn new(xmin: f32, xmax: f32, ymin: f32, ymax: f32) -> Self {
        let mut ret = Self {
            vert: Default::default(),
            init_vert: 0,
        };

        // Filled counterclockwise
        // (3)xmin,ymax----------(2)xmax,ymax
        //  |                     |
        //  |                     |
        //(0)xmin,ymin-----------(1)xmax,ymin
        
        ret.vert.push(Vertex{pos: Vector2::new(xmin, ymin), ccw: 1, cw: 3});
        ret.vert.push(Vertex{pos: Vector2::new(xmax, ymin), ccw: 2, cw: 0});
        ret.vert.push(Vertex{pos: Vector2::new(xmax, ymax), ccw: 3, cw: 1});
        ret.vert.push(Vertex{pos: Vector2::new(xmin, ymax), ccw: 0, cw: 2});

        ret
    }

    #[inline(always)]
    fn dist(&self, i: usize, line: &CuttingLine) -> f32 {
        (line.pos.dot(&self.vert[i].pos)-line.r*line.r) / line.r
    }

    #[inline(always)]
    fn next(&self, i: usize) -> usize {
        self.vert[i].ccw
    }

    #[inline(always)]
    fn pos(&self, i: usize) -> &Vector2 {
        &self.vert[i].pos
    }

    #[inline(always)]
    fn vertex(&mut self, i: usize) -> &mut Vertex {
        &mut self.vert[i]
    }

    pub fn add_point(&mut self,point: &Vector2) -> bool {
        let line = CuttingLine::new(point);
        let cut1;
        let cut2;

        // Find first point which is inside
        // It is guaranteed that there is one, no need to check return to zero
        let mut cur_i = self.init_vert; // Start from init
        let mut cur_d = self.dist(cur_i, &line);
        while cur_d >= TOL {
            //self.to_delete.push(cur_i);
            cur_i = self.next(cur_i);
            cur_d = self.dist(cur_i, &line)
        }

        //println!("First inside point: {cur_i}");

        // Now we are at first inside point. Update init point
        self.init_vert = cur_i;

        // Search for first outside point
        loop {
            let next_i = self.next(cur_i);
            if next_i==self.init_vert {
                // All points are inside, nothing to do
                return false;
            }
            let next_d = self.dist(next_i, &line);
            if next_d >= TOL {
                // next_i is the first outside point
                // cur_i is the last inside point
                cut1 = CutEdge {
                    in_ind: cur_i,
                    in_dist: cur_d,
                    out_ind: next_i,
                    out_dist: next_d
                };
                // Advance to first outside point
                cur_i = next_i;
                cur_d = next_d;
                break;
            }
            cur_i = next_i;
            cur_d = next_d;
        }

       // println!("First outside point: {cur_i}");

        // Now cur_i is the first outside point
        // Search again for the next inside point
        loop {
            let next_i = self.next(cur_i);
            let next_d = self.dist(next_i, &line);
            if next_d < TOL {
                // next_i is the first inside point
                // cur_i is the last outside point
                cut2 = CutEdge {
                    out_ind: cur_i,
                    out_dist: cur_d,
                    in_ind: next_i,
                    in_dist: next_d
                };
                //cur_i = next_i;
                //cur_d = next_d;
                break;
            }
            cur_i = next_i;
            cur_d = next_d;
        }
        
        //println!("Second inside point: {cur_i}");
        //println!("Cut 1: {cut1:?}");
        //println!("Cut 1: {cut2:?}");

        // At this point both cuts are found.
        
        // Cut #1:
        // For first cut we can always reuse the outer point
        // But we need to defer modifying the point coordinates
        // Until second cut is done because they may share the same outer point
        let cut1_point = self.vertex_pos_from_cut(&cut1);
        // Connection to cut2 will be set later
        
        // Cut #2:
        let p = self.vertex_pos_from_cut(&cut2);
        // Check if we can reuse the outer point (should not be the same point as for cut1)
        if cut1.out_ind != cut2.out_ind {
            let v = self.vertex(cut2.out_ind);
            v.pos = p;
            v.cw = cut1.out_ind; // Link to vertex from first cut
            // Set missing ccw connection cut1->cut2
            self.vertex(cut1.out_ind).ccw = cut2.out_ind;
        } else {
            // There is only one outer point, so we need to add new vertex
            self.vert.push(
                Vertex {
                    pos: p,
                    cw: cut1.out_ind, // Link to vertes from first cut
                    ccw: cut2.in_ind, // Reset connection to next inner point
                }
            );
            // Set missing ccw connection cut1->cut2
            self.vertex(cut1.out_ind).ccw = self.vert.len()-1; // New added point
        }

        // Defered update of outer point in cut1
        self.vertex(cut1.out_ind).pos = cut1_point;
        
        true
    }

    #[inline(always)]
    fn vertex_pos_from_cut(&self, cut: &CutEdge) -> Vector2 {
        let frac = cut.out_dist / (cut.in_dist.abs()+cut.out_dist);
        // New point
        //self.pos(cut.out_ind)+frac*(self.pos(cut.in_ind)-self.pos(cut.out_ind))
        (1.0-frac)*self.pos(cut.out_ind)+frac*self.pos(cut.in_ind)
    }

    pub fn write_to_file(&self, fname: &str) {
        let mut s = String::new();
        let mut i = self.init_vert;
        loop {
            s.push_str(&format!("{} {}\n",self.vert[i].pos[0],self.vert[i].pos[1]));
            i = self.vert[i].ccw;
            if i == self.init_vert {break}
        }
        let mut out = File::create(fname).unwrap();
        write!(out,"{}",s).unwrap();
    }

}

#[cfg(test)]
mod tests {
    use super::{Vector2, VoronoiCell};

    #[test]
    fn cell_test() {
        let mut c = VoronoiCell::new(-10.0, 10.0, -10.0, 10.0);
        c.add_point(&Vector2::new(20.0,-17.5));
        c.add_point(&Vector2::new(10.0,0.0));
        c.add_point(&Vector2::new(10.0,10.0));
        c.add_point(&Vector2::new(-5.0*2.0,-2.0*2.0));
        c.add_point(&Vector2::new(-5.0*2.0,5.0*2.0));
        c.write_to_file("tests/out.dat");
    }
}