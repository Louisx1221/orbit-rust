//@file       : math.rs
//@autor      : github.com/louisx1221
//@date       : 2021/11/22

//use std::io;
//use rand::Rng;
mod math;
use crate::math::*;
mod orbit;
use crate::orbit::*;

fn main() {
    let mut orb = Orb {
        r: [731.8353302913914, -5059.425707600824, 4610.677207629951],
        v: [-0.5747084918184473, -5.158344977093602, -5.564664141587914],
        coe: [0.; 6],
        jd: gre2julian([2021., 7., 1., 06., 0., 0.]),
    };
    orb.rv2coe();
    let rv = [orb.r[0], orb.r[1], orb.r[2], orb.v[0], orb.v[1], orb.v[2]];
    orb.orb_prop(86400.);
    orb.coe2rv();
    println!("orb is {:?}", orb);
    let gd = julian2gre(orb.jd);
    println!("gd is {:?}", gd);
    let x = rkf78(orb_prop_j2, rv.to_vec(), 0., 86400., 1.);
    println!("x is {:?}", x);
}