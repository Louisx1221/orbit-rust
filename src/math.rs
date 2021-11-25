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

pub fn unit(vec: [f64; 3]) -> [f64; 3] {
    let vec_mag = norm(vec);
    let mut out = vec;
    for element in out.iter_mut() {
        *element /= vec_mag;
    }
    out
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

pub fn rkf78(eq: fn(Vec<f64>) -> Vec<f64>, x0: Vec<f64>, ti: f64, tf: f64, h: f64) -> Vec<f64> { 
    // Runge-Kutta-Fehlberg
    // Fehlberg's 7th and 8th Order Embedded Method

    // define integration coefficients
    // set initialization indicator
    static CH: [f64; 13] = [0., 0., 0., 0., 0., 34. / 105., 9. / 35., 9. / 35., 9. / 280., 9. / 280., 0., 41. / 840., 41. / 840.];
	static ALPHA: [f64; 13] = [0., 2. / 27., 1. / 9., 1. / 6., 5. / 12., 0.5, 5. / 6., 1. / 6., 2. / 3., 1. / 3., 1., 0., 1.];
	static BETA: [[f64; 12]; 13] = [[0.,             0.,       0.,         0.,           0.,            0.,          0.,            0.,        0.,         0.,        0., 0., ],
                                    [2. / 27.,       0.,       0.,         0.,           0.,            0.,          0.,            0.,        0.,         0.,        0., 0., ],
                                    [1. / 36.,       1. / 12., 0.,         0.,           0.,            0.,          0.,            0.,        0.,         0.,        0., 0., ],
                                    [1. / 24.,       0.,       1. / 8.,    0.,           0.,            0.,          0.,            0.,        0.,         0.,        0., 0., ],
                                    [5. / 12.,       0.,       -25. / 16., 25. / 16.,    0.,            0.,          0.,            0.,        0.,         0.,        0., 0., ],
                                    [0.05,           0.,       0.,         0.25,         0.2,           0.,          0.,            0.,        0.,         0.,        0., 0., ],
                                    [-25. / 108.,    0.,       0.,         125. / 108.,  -65. / 27.,    125. / 54.,  0.,            0.,        0.,         0.,        0., 0., ],
                                    [31. / 300.,     0.,       0.,         0.,           61. / 225.,    -2. / 9.,    13. / 900.,    0.,        0.,         0.,        0., 0., ],
                                    [2.,             0.,       0.,         -53. / 6.,    704. / 45.,    -107. / 9.,  67. / 90.,     3.,        0.,         0.,        0., 0., ],
                                    [-91. / 108.,    0.,       0.,         23. / 108.,   -976. / 135.,  311. / 54.,  -19. / 60.,    17. / 6.,  -1. / 12.,  0.,        0., 0., ],
                                    [2383. / 4100.,  0.,       0.,         -341. / 164., 4496. / 1025., -301. / 82., 2133. / 4100., 45. / 82., 45. / 164., 18. / 41., 0., 0., ],
                                    [3. / 205.,      0.,       0.,         0.,           0.,            -6. / 41.,   -3. / 205.,    -3. / 41., 3. / 41.,   6. / 41.,  0., 0., ],
                                    [-1777. / 4100., 0.,       0.,         -341. / 164., 4496. / 1025., -289. / 82., 2193. / 4100., 51. / 82., 33. / 164., 12. / 41., 0., 1., ]];

    // let mut x = Vec::new();
    // for i in 0..x.len() {
    //     x.push(x0[i]);
    // }
    let mut x = x0.clone();
    let mut xwrk;
    let mut xdot;
    let neq = x.len();
    let mut f: Vec<[f64; 13]> = Vec::new();
    for _i in 0..neq {
        f.push([0.; 13]);
    }

    let mut t0 = ti;
    let tetol = 1e-12;
    let mut twrk;
    let mut xerr;
    let mut ter;
    let mut tol;
    let mut tconst;
    let mut ab;
    let mut ba;

    //compute integration "direction"
	let sdt = (tf - ti).signum();
	let mut dt = h.abs() * sdt;

    loop {
        // load "working" time and integration vector
        twrk = t0;
        xwrk = x.clone();
        // check for last dt
        if dt.abs() > (tf - t0).abs() {
			dt = tf - t0;
        }
		if (t0 - tf).abs() < 0.00000001 {
			break
		}
		// evaluate equations of motion
		//xdot = feval(deq, ti, x);
		xdot = eq(x.clone());
		for i in 0..neq {
			f[i][0] = xdot[i];
        }
		// compute solution
		for k in 1..13
		{
			for i in 0..neq {
                x[i] = 0.;
				for j in 0..k {
					x[i] += BETA[k][j] * f[i][j];
                }
                x[i] *= dt;
				x[i] += xwrk[i];
			}
			t0 = twrk + ALPHA[k] * dt;
			//xdot = feval(deq, ti, x);
			xdot = eq(x.clone());
			for j in 0..neq {
				f[j][k] = xdot[j];
            }
		}
		for i in 0..neq {
			x[i] = 0.;
			for l in 0..13 {
				x[i] += CH[l] * f[i][l];
            }
            x[i] *= dt;
            x[i] += xwrk[i];
		}

		// truncation error calculations
		xerr = tetol;
		for i in 0..neq {
			ter = ((f[i][0] + f[i][10] - f[i][11] - f[i][12]) * CH[11] * dt).abs();
			tol = x[i].abs() * tetol + tetol;
			tconst = ter / tol;
			if tconst > xerr {
				xerr = tconst;
            }
		}

		// compute new step size
		//dt = 0.8 * dt * (1 / xerr) ^ (1 / 8);
		ab = 1. / xerr;
		ba = 1. / 8.;
        dt *= ab.powf(ba);
        dt *= 0.8;
		if xerr > 1. {
			// reject current step
			t0 = twrk;
			x = xwrk.clone();
		}

    }

    x
}