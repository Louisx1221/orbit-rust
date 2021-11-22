//use std::io;
//use rand::Rng;
mod math;
// pub use crate::math::*;
mod orbit;
use crate::orbit::*;

fn main() {
    let mut orb = Orb {
        r: [-4.7586e6, 1.7409e6, 3.8969e6],
        v: [-4.9457e3, -6.7309e3, -0.9616e3],
        coe: [6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
        jd: 12.0,
    };
    //println!("r = {} {} {}", orb.r[0], orb.r[1], orb.r[2]);
    orb.rv2coe();
    println!("orb is {:#?}", orb);

    let ea = solve_kepler(2.0, 0.2);
    println!("ea = {}", ea);
}