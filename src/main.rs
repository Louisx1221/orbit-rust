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

    // let ri_mag = 7000.;
    // let rf_mag = 8000.;
    // let mut ri = [ri_mag, 1., 1.];
    // let mut rf = [-rf_mag, -1., -1.];
    // let vt = hohnmann_trans(ri_mag, rf_mag);
    // println!("vt is {:?}", vt);
    // let mut sol = lambert_izzo2015(ri, rf, vt.1, 1., 0);
    // println!("sol is {:?}", sol);
    // let rv = [ri[0], ri[1], ri[2], sol.0[0][0] as f64, sol.0[0][1] as f64, sol.0[0][2] as f64];
    // let x = rkf78(orb_prop_j2, rv.to_vec(), 0., vt.1, 1.);
    // println!("x is {:?}", x);
    // vt is ([0.24747702649420322, 0.23934746429774822, 0.48682449079195145], 3232.0114915803797)
    // sol is ([[-0.0015221994949434652, 5.5108577678602595, 5.5108577678602595]], [[0.0012455304284963592, -4.822000581461361, -4.822000581461361]], 1)
    // x is [-7985.436604402226, -18.283477773160186, -39.333129921757134, 0.03682691415518792, -4.8307105151758485, -4.830094801381817]

    let ri =[5000., 10000., 2100.];
    let rf = [-14600., 2500., 7000.];
    let tof = 3600.;
    let sol = lambert_izzo2015(ri, rf, tof, 1., 0);
    println!("sol is {:?}", sol);
    // v1 (km/s) = [-5.99249 1.92536 3.24564]
    // v2 (km/s) = [-3.31246 -4.19662 -0.385288]
}