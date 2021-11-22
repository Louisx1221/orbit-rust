//@file       : main.rs
//@autor      : github.com/louisx1221
//@date       : 2021/11/22

pub use core::f64::consts::PI;

pub const TWO_PI: f64 = 2. * PI;

pub fn fmod(val_in: f64, val_mod: f64) -> f64 {
    let mut val_out = val_in;
    while val_out > val_mod {
        val_out -= val_mod;
    }
    while val_out < 0. {
        val_out += val_mod;
    }
    val_out
}

pub fn norm(vec: [f64; 3]) -> f64 {
    (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]).sqrt()
}

pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let mut c: [f64; 3] = [0.0; 3];
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    c
}

pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let mut c: f64 = 0.;
    for i in 0..3 {
        c += a[i] * b[i];
    }
    c
}

pub fn atan2(y: f64, x: f64) -> f64 {
    let mut out = (y / x).atan();
    out = if x < 0. {
        if y > 0. {
            out + PI
        } else {
            out - PI
        }
    } else {
        out
    };
    out
}
