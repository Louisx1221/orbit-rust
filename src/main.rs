//@file       : main.rs
//@autor      : github.com/louisx1221
//@date       : 2021/11/22

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
        coe: [6., 7., 8., 9., 10., 11.],
        jd: 12.,
    };
    orb.rv2coe();
    orb.orb_prop(86400.);
    orb.coe2rv();
    println!("orb is {:#?}", orb);
}