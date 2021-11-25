//@file       : math.rs
//@autor      : github.com/louisx1221
//@date       : 2021/11/22

//use std::io;
//use rand::Rng;
mod math;
use crate::math::*;
mod orbit;
use crate::orbit::*;
mod lambert;
use crate::lambert::*;

fn main() {
    // let mut orb = Orb {
    //     r: [731.8353302913914, -5059.425707600824, 4610.677207629951],
    //     v: [-0.5747084918184473, -5.158344977093602, -5.564664141587914],
    //     coe: [0.; 6],
    //     jd: gre2julian([2021., 7., 1., 06., 0., 0.]),
    // };
    // orb.rv2coe();
    // let rv = [orb.r[0], orb.r[1], orb.r[2], orb.v[0], orb.v[1], orb.v[2]];
    // orb.orb_prop(86400.);
    // orb.coe2rv();
    // println!("orb is {:?}", orb);
    // let gd = julian2gre(orb.jd);
    // println!("gd is {:?}", gd);
    // let x = rkf78(orb_prop_j2, rv.to_vec(), 0., 86400., 1.);
    // println!("x is {:?}", x);

    let ri_mag = 7000.;
    let rf_mag = 8000.;
    let ri = [ri_mag, 1., 1.];
    let rf = [-rf_mag, -1., -1.];
    let vt = hohnmann_trans(ri_mag, rf_mag);
    println!("vt is {:?}", vt);
    let sol = lambert_izzo2015(ri, rf, vt.1, 1., 0);
    println!("sol is {:?}", sol);
    let rv = [ri[0], ri[1], ri[2], sol.0[0][0] as f64, sol.0[0][1] as f64, sol.0[0][2] as f64];
    let x = rkf78(orb_prop_j2, rv.to_vec(), 0., vt.1, 1.);
    println!("x is {:?}", x);
}