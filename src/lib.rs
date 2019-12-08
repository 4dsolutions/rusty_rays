pub mod vlib {

    pub struct Vivm {
        pub a: f64,
        pub b: f64,
        pub c: f64,
        pub d: f64,
    }

    pub struct Vxyz {
        pub x: f64,
        pub y: f64,
        pub z: f64,
    }

    impl Vxyz {
        pub fn to_ivm(&self) -> Vivm {
            let (x,y,z)    =  (self.x, self.y, self.z);
            let k:f64        =  2.0/2.0_f64.sqrt();

            let mut xge0: f64 = 0.0;
            let mut xlt0: f64 = 0.0;
            let mut yge0: f64 = 0.0;
            let mut ylt0: f64 = 0.0;
            let mut zge0: f64 = 0.0;
            let mut zlt0: f64 = 0.0;

            if x >= 0.0 {xge0 = 1.0} else {xlt0 = 1.0};
            if y >= 0.0 {yge0 = 1.0} else {ylt0 = 1.0};
            if z >= 0.0 {zge0 = 1.0} else {zlt0 = 1.0};
  
            let a = k * (xge0 * ( x) + yge0 * ( y) + zge0 * ( z));
            let b = k * (xlt0 * (-x) + ylt0 * (-y) + zge0 * ( z));
            let c = k * (xlt0 * (-x) + yge0 * ( y) + zlt0 * (-z));
            let d = k * (xge0 * ( x) + ylt0 * (-y) + zlt0 * (-z));

            Vivm{ 
                a: a,
                b: b,
                c: c,
                d: d,
            }
        }           
    }

    impl Vivm {

        pub fn to_xyz(&self) -> Vxyz {
            let (a,b,c,d)    =  (self.a, self.b, self.c, self.d);
            let k:f64        =  0.5/2.0_f64.sqrt();
            let xyz:[f64; 3] = [k * (a - b - c + d),
                                k * (a - b + c - d),
                                k * (a + b - c - d)];
            Vxyz{
                x: xyz[0],
                y: xyz[1],
                z: xyz[2]
            }           
        }
    }

    pub trait Ray {
        fn new(coords: &[f64]) -> Self;
        fn add(&self, other: &Self) -> Self;
        fn neg(&self) -> Self;
        fn sub(&self, other: &Self) -> Self;
        fn mul(&self, scalar: f64) -> Self;
        fn dot(&self, other: &Self) -> f64;
        fn length(&self) -> f64;
    }
 
    impl Ray for Vivm {
        fn new(coords: &[f64]) -> Vivm {
            let mut the_min = coords[0];
            for numb in coords[1..].iter() {
                the_min = the_min.min(*numb)
            }
           
            let new_coords:[f64;4] = [coords[0] - the_min, 
                                      coords[1] - the_min,
                                      coords[2] - the_min,
                                      coords[3] - the_min];

            Vivm{
            a: new_coords[0],
            b: new_coords[1],
            c: new_coords[2],
            d: new_coords[3]
            }
        }

        fn add(&self, other: &Vivm) -> Vivm {
            Vivm::new(&[
                self.a + other.a,
                self.b + other.b,
                self.c + other.c,
                self.d + other.d])
        }

        fn neg(&self) -> Vivm {
            Vivm::new(&[
                -self.a,
                -self.b,
                -self.c,
                -self.d])

        }

        fn sub(&self, other: &Vivm) -> Vivm {
            Vivm::new(&[
                self.a - other.a,
                self.b - other.b,
                self.c - other.c,
                self.d - other.d])
        }

        fn mul(&self, scalar: f64) -> Vivm {
            Vivm{
                a: self.a * scalar,
                b: self.b * scalar,
                c: self.c * scalar,
                d: self.d * scalar}
        }

        fn dot(&self, other: &Vivm) -> f64 {
            // s1 = a.dot(b)/(a.length() * b.length())
            // degrees(acos(s1))
            // 109.47122063449069
            let avg = (self.a + self.b + self.c + self.d)/4.0;
            let (sa, sb, sc, sd) = (self.a - avg, self.b - avg, self.c - avg, self.d - avg);
            let avg = (other.a + other.b + other.c + other.d)/4.0;
            let (oa, ob, oc, od) = (other.a - avg, other.b - avg, other.c - avg, other.d - avg);
            return 0.5 * (sa*oa + sb*ob + sc*oc + sd*od);
        }

        fn length(&self) -> f64 {
            // Return this vector's length"""
            return self.dot(&self).sqrt();
        }

    }

    impl Ray for Vxyz {

        fn new(coords: &[f64]) -> Vxyz {
            Vxyz{
            x: coords[0],
            y: coords[1],
            z: coords[2]
            }
        }

        fn add(&self, other: &Vxyz) -> Vxyz {
            Vxyz{
                x: self.x + other.x,
                y: self.y + other.y,
                z: self.z + other.z}
        }

        fn neg(&self) -> Vxyz {
            Vxyz{
                x: -self.x,
                y: -self.y,
                z: -self.z}
        }

        fn sub(&self, other: &Vxyz) -> Vxyz {
            Vxyz{
                x: self.x - other.x,
                y: self.y - other.y,
                z: self.z - other.z}
        }

        fn mul(&self, scalar: f64) -> Vxyz {
            Vxyz{
                x: self.x * scalar,
                y: self.y * scalar,
                z: self.z * scalar
            }

        }

        fn dot(&self, other: &Vxyz) -> f64 {
            let (sx,sy,sz) = (self.x, self.y, self.z);
            let (ox,oy,oz) = (other.x, other.y, other.z);
            return sx * ox + sy * oy + sz * oz;
        }

        fn length(&self) -> f64 {
            // Return this vector's length"""
            return self.dot(&self).sqrt();
        }

    }
}

pub mod tetrahedron {

    pub struct TetEdges {
        pub ab: f64,
        pub ac: f64,
        pub ad: f64,
        pub bc: f64,
        pub cd: f64,
        pub db: f64,
    }

    impl TetEdges {
        pub fn volume(&self) -> f64 {
            let a = self.ab;
            let b = self.ac;
            let c = self.ad;
            let d = self.bc;
            let e = self.cd;
            let f = self.db;
            let a2 = a*a;
            let b2 = b*b;
            let c2 = c*c;
            let d2 = d*d;
            let e2 = e*e;
            let f2 = f*f;

            let mut sumval = f2 * a2 * b2;
            sumval +=  d2 * a2 * c2;
            sumval +=  a2 * b2 * e2;
            sumval +=  c2 * b2 * d2;
            sumval +=  e2 * c2 * a2;
            sumval +=  f2 * c2 * b2;
            sumval +=  e2 * d2 * a2;
            sumval +=  b2 * d2 * f2;
            sumval +=  b2 * e2 * f2;
            sumval +=  d2 * e2 * c2;
            sumval +=  a2 * f2 * e2;
            sumval +=  d2 * f2 * c2; 
            let addopen = sumval;

            let mut sumval = a2 * b2 * d2;
            sumval +=  d2 * e2 * f2;
            sumval +=  b2 * c2 * e2;
            sumval +=  a2 * c2 * f2;
            let addclosed = sumval;   
            
            let mut sumval =  a2 * e2 * (a2 + e2);
            sumval += b2 * f2 * (b2 + f2);
            sumval += c2 * d2 * (c2 + d2);
            let addopposite = sumval;

            ((addopen - addclosed - addopposite)/2.0).sqrt()
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use vlib::Ray;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn cube() {

        let unit_tet = tetrahedron::TetEdges {
            ab: 1.0,
            ac: 1.0,
            ad: 1.0,
            bc: 1.0,
            cd: 1.0,
            db: 1.0,
        };

        let corner_tet = tetrahedron::TetEdges {
            ab: 2_f64.sqrt()/2.0,
            ac: 2_f64.sqrt()/2.0,
            ad: 2_f64.sqrt()/2.0,
            bc: 1.0,
            cd: 1.0,
            db: 1.0,
        };
        assert_eq!(4.0  * corner_tet.volume() +
                                unit_tet.volume(), 3.0)
    }

    #[test]
    fn add_xyz(){
        let x1 = vlib::Vxyz::new(&[1.0_f64, 1.0, 1.0]);
        let x2 = vlib::Vxyz::new(&[1.0_f64, 1.0, 1.0]);
        let r = x1.add(&x2); 
        assert_eq!((2.0_f64, 2.0, 2.0), (r.x, r.y, r.z));
    }
}