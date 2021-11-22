pub use core::f32::consts::PI;

pub const TWO_PI: f32 = 2.0 * PI;

pub fn fmod(val_in: f32, val_mod: f32) -> f32 {
    let mut val_out = val_in;
    while val_out > val_mod {
        val_out -= val_mod;
    }
    while val_out < 0.0 {
        val_out += val_mod;
    }
    val_out
}

pub fn norm(vec: [f32; 3]) -> f32 {
    (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]).sqrt()
}

pub fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    let mut c: [f32; 3] = [0.0; 3];
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    c
}

pub fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    let mut c: f32 = 0.0;
    for i in 0..3 {
        c += a[i] * b[i];
    }
    c
}

pub fn atan2(y: f32, x: f32) -> f32 {
    let mut out = (y / x).atan();
    out = if x < 0.0 {
        if y > 0.0 {
            out + PI
        } else {
            out - PI
        }
    } else {
        out
    };
    out
}
