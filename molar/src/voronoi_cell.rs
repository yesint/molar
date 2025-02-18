use std::fs::File;
use std::io::Write;

pub type Vector2f = nalgebra::Vector2<f32>;
const TOL: f32 = 1e-10;

#[derive(Debug)]
pub struct Vertex {
    pos: Vector2f,
    // Index of neighbour in ccw direction
    ccw_neib: usize, 
    // ID of point that created the ccw edge
    // Negative numbers are reserved for the walls
    edge_id: i32, 
}

struct CuttingLine {
    pos: Vector2f,
    r2: f32,
}

impl CuttingLine {
    fn new(point: &Vector2f) -> Self {
        let pos = 0.5*point;
        let r2 = pos.norm_squared();
        Self {pos, r2}
    }
}

pub struct VoronoiCell {
    // Vertices
    vert: Vec<Vertex>,
    init_vert: usize,
}

#[derive(Default,Debug)]
struct EdgeCut {
    in_ind: usize,
    in_dist: f32,
    out_ind: usize,
    out_dist: f32,
}

impl VoronoiCell {
    pub fn new(xmin: f32, xmax: f32, ymin: f32, ymax: f32) -> Self {
        let mut ret = Self {
            vert: Vec::with_capacity(8),
            init_vert: 0,
        };

        // Filled counterclockwise
        // (3)xmin,ymax----------(2)xmax,ymax
        //  |                     |
        //  |                     |
        //(0)xmin,ymin-----------(1)xmax,ymin
        
        ret.vert.push(Vertex{pos: Vector2f::new(xmin, ymin), ccw_neib: 1, edge_id: -1});
        ret.vert.push(Vertex{pos: Vector2f::new(xmax, ymin), ccw_neib: 2, edge_id: -2});
        ret.vert.push(Vertex{pos: Vector2f::new(xmax, ymax), ccw_neib: 3, edge_id: -3});
        ret.vert.push(Vertex{pos: Vector2f::new(xmin, ymax), ccw_neib: 0, edge_id: -4});

        ret
    }

    #[inline(always)]
    fn dist(&self, i: usize, line: &CuttingLine) -> f32 {
        line.pos.dot(&self.pos(i))-line.r2
    }

    #[inline(always)]
    fn next(&self, i: usize) -> usize {
        unsafe{ self.vert.get_unchecked(i).ccw_neib }
    }

    #[inline(always)]
    fn pos(&self, i: usize) -> &Vector2f {
        unsafe{ &self.vert.get_unchecked(i).pos }
    }

    #[inline(always)]
    fn vertex(&self, i: usize) -> &Vertex {
        unsafe{ self.vert.get_unchecked(i) }
    }

    #[inline(always)]
    fn vertex_mut(&mut self, i: usize) -> &mut Vertex {
        unsafe{ self.vert.get_unchecked_mut(i) }
    }

    pub fn add_point(&mut self,point: &Vector2f, id: usize) -> bool {
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
                cut1 = EdgeCut {
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

        // Now cur_i is the first outside point
        // Search again for the next inside point
        loop {
            let next_i = self.next(cur_i);
            let next_d = self.dist(next_i, &line);
            if next_d < TOL {
                // next_i is the first inside point
                // cur_i is the last outside point
                cut2 = EdgeCut {
                    out_ind: cur_i,
                    out_dist: cur_d,
                    in_ind: next_i,
                    in_dist: next_d
                };
                break;
            }
            cur_i = next_i;
            cur_d = next_d;
        }
        
        // At this point both cuts are found.
        
        // Cut #2:
        let p = self.vertex_pos_from_cut(&cut2);
        // Check if we can reuse the outer point 
        // (should not be the same point as for cut1)
        if cut1.out_ind != cut2.out_ind {
            self.vertex_mut(cut2.out_ind).pos = p;
            // Set ccw connection cut1->cut2
            // For first cut we always reuse its outer point
            self.vertex_mut(cut1.out_ind).ccw_neib = cut2.out_ind;
        } else {
            // There is only one outer point, so we need to add new vertex
            let old_id = self.vertex(cut2.out_ind).edge_id;
            self.vert.push(
                Vertex {
                    pos: p,
                    ccw_neib: cut2.in_ind, // Reset connection to next inner point
                    edge_id: old_id // Reset ID
                }
            );
            // Set ccw connection cut1->cut2
            // For first cut we always reuse its outer point
            self.vertex_mut(cut1.out_ind).ccw_neib = self.vert.len()-1; // New added point
        }

        // Cut #1:
        // For first we always reuse the outer point
        // Second cut is already done so no problem if 
        // they share the same outer point
        self.vertex_mut(cut1.out_ind).pos = self.vertex_pos_from_cut(&cut1);
        self.vertex_mut(cut1.out_ind).edge_id = id as i32;
        
        true
    }

    #[inline(always)]
    fn vertex_pos_from_cut(&self, cut: &EdgeCut) -> Vector2f {
        let frac = cut.out_dist / (cut.in_dist.abs()+cut.out_dist);
        (1.0-frac)*self.pos(cut.out_ind)+frac*self.pos(cut.in_ind)
    }

    pub fn iter_vertex(&self) -> VoronoiCellVertexIter {
        VoronoiCellVertexIter{
            cell: self,
            cur: self.init_vert,
            stop: false
        }
    }

    pub fn write_to_file(&self, fname: &str) {
        let mut s = String::new();
        self.for_each_vert(
            |v| s.push_str(&format!("{} {}\n",v.pos[0],v.pos[1]))
        );
        let mut out = File::create(fname).unwrap();
        write!(out,"{}",s).unwrap();
    }

    pub fn area(&self) -> f32 {
        let mut a = 0.0;
        let mut i = self.init_vert;
        loop {
            let v1 = self.vertex(i);
            i = v1.ccw_neib;
            let v2 = self.vertex(i);
            a += 0.5*v1.pos.cross(&v2.pos).norm();
            if i == self.init_vert {break}
        }
        a
    }

    pub fn for_each_vert(&self, mut func: impl FnMut(&Vertex)) {
        let mut i = self.init_vert;
        loop {
            let v = self.vertex(i);
            func(v);
            i = v.ccw_neib;
            if i == self.init_vert {break}
        }
    }

}

pub struct VoronoiCellVertexIter<'a> {
    cell: &'a VoronoiCell,
    cur: usize,
    stop: bool,
}

impl<'a> Iterator for VoronoiCellVertexIter<'a> {
    type Item = &'a Vertex;
    fn next(&mut self) -> Option<Self::Item> {
        if self.stop {
            None
        } else {
            let ret = &self.cell.vert[self.cur];
            self.cur = ret.ccw_neib;
            if self.cur == self.cell.init_vert {
                self.stop = true;
            }
            Some(ret)
        }
    }
}


#[cfg(test)]
mod tests {
    use super::{Vector2f, VoronoiCell};

    #[test]
    fn cell_test() {
        let mut c = VoronoiCell::new(-10.0, 10.0, -10.0, 10.0);
        c.add_point(&Vector2f::new(20.0,-17.5), 1);
        c.add_point(&Vector2f::new(10.0,0.0), 2);
        c.add_point(&Vector2f::new(10.0,10.0), 3);
        c.add_point(&Vector2f::new(-5.0*2.0,-2.0*2.0), 4);
        c.add_point(&Vector2f::new(-5.0*2.0,5.0*2.0), 5);
        c.write_to_file("target/out.dat");

        for v in c.iter_vertex() {
            println!("{v:?}")
        }
    }
}