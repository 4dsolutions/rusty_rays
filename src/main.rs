use rays::vlib;
use rays::vlib::Ray;

fn main() {
    println!("Hello, world");
    let coords = [1.0_f64, 1.0, 1.0];
    let qa = vlib::Vxyz::new(&coords);
    println!("{} {} {}", qa.x, qa.y, qa.z);
    let coords = [2.0_f64, 1.0, 1.0, 0.0];
    let qb = vlib::Vivm::new(&coords);
    println!("{} {} {} {}", qb.a, qb.b, qb.c, qb.d);

    let x1 = vlib::Vxyz::new(&[1.0_f64, 1.0, 1.0]);
    let x2 = vlib::Vxyz::new(&[1.0_f64, 1.0, 1.0]);
    let r = x1.add(&x2); 
    println!("{} {} {}", r.x, r.y, r.z);

    let negqb = qb.neg();
    println!("{} {} {} {}", negqb.a, negqb.b, negqb.c, negqb.d);
   
    let xyz = negqb.to_xyz();
    println!("{} {} {}", xyz.x, xyz.y, xyz.z); 
    
    let ivm = xyz.to_ivm();
    println!("{} {} {} {}", ivm.a, ivm.b, ivm.c, ivm.d);

    println!("Length of ivm = {}", ivm.length());
    println!("Length of xyz = {}", xyz.length());
}