//@file       : lambert.rs
//@autor      : github.com/louisx1221
//@date       : 2021/11/25

use crate::math::*;
use crate::orbit::*;

pub fn lambert_izzo2015(r1: [f64; 3], r2: [f64; 3], tof: f64, cw: f64, multi_revs: u32)-> (Vec<[f64; 3]>, Vec<[f64; 3]>, u32) {
    // solve Lambert's orbital two point boundary value problem (2015 Dario Izzo)

    // input
    // r1: starting position (x1,y1,z1)
    // r2: final position (x2,y2,z2)
    // tof: time of flight
    // cbmu: gravitational parameter 
    // cw: 1 for retrograde motion (clockwise), -1 if counter-clock wise
    // max_revs: Maximum number of multirevs to be computed 

    // output
    // nsol: number of solutions
    // m_v1: matrix of state vector solutions of initial velocity vector of the transfer orbit
    // m_v2: matrix of state vector solutions of final velocity vector of the transfer orbit

    // m_v1 or m_v2(sn, 1) : velocity vector x component
    // m_v1 or m_v2(sn, 2) : velocity vector y component
    // m_v1 or m_v2(sn, 3) : velocity vector z component

    // where sn is the solution number
    // toc
    // tic

    // global m_lambda

    // 0 - Sanity checks
    if tof <= 0. {
        println!("Time of flight is negative!");
    }

    // 1 - Getting lambda and T
    let m_c = ((r2[0] - r1[0]) * (r2[0] - r1[0]) + (r2[1] - r1[1]) * (r2[1] - r1[1]) + (r2[2] - r1[2]) * (r2[2] - r1[2])).sqrt();
    let r1_mag = norm(r1);
    let r2_mag = norm(r2);
    let m_s = (m_c + r1_mag + r2_mag) / 2.;

    let ir1 = unit(r1);
    let ir2 = unit(r2);
    
    let mut ih = cross(ir1, ir2);
    ih = unit(ih);

    let lambda2 = 1. - m_c / m_s;
    let mut m_lambda = lambda2.sqrt();

    let mut it1: [f64; 3];
    let mut it2: [f64; 3];
    if ih[2] * cw < 0. {
        // Transfer angle is larger than 180 degrees as seen from abive the z axis
        m_lambda = -m_lambda;
        it1 = cross(ir1, ih);
        it2 = cross(ir2, ih);
    } else { 
        it1 = cross(ih, ir1);
        it2 = cross(ih, ir2);
    }
    it1 = unit(it1);
    it2 = unit(it2);

    let lambda3 = m_lambda * lambda2;
    let k_p = (2. * MU / m_s / m_s / m_s).sqrt();
    let t = k_p * tof;

    // 2 - We now have lambda, T and we will find all x
    // 2.1 - Let us first detect the maximum number of revolutions for which there exists a solution
    let mut m_nmax = (t / PI).floor() as u32;
    let t00 = m_lambda.acos() + m_lambda * (1. - lambda2).sqrt();
    let t0 = t00 + m_nmax as f64 * PI;
    let t1 = 2. / 3. * (1. - lambda3);

    let mut it;
    let mut t_min;
    let mut x_old;
    let mut x_new;
    let mut err; 
    let mut dt;
    let mut ddt;
    let mut dddt;
    let mut dt3;

    if m_nmax > 0 {
        if t < t0 {
            // We use Halley iterations to find xM and TM
            it = 0;
            //  err = 1.;
            t_min = t0;
            x_old = 0.;
            x_new = 0.;
            loop {
                dt3 = dtdx_izzo(m_lambda, x_old, t_min);
                dt = dt3.0;
                ddt = dt3.1;
                dddt = dt3.2;

                if dt != 0. {
                    x_new = x_old - dt * ddt / (ddt * ddt - dt * dddt / 2.);
                }
                err = (x_old - x_new).abs();
                if (err < 1e-13) || (it > 50) { 
                    break
                }
                t_min = x2tof_izzo(m_lambda, x_new, m_nmax);
                x_old = x_new;
                it += 1;
            }
            if t_min > t {
                m_nmax -= 1;
            }
        }
    }
    // We exit this if clause with Mmax being the maximum number of revolutions
    // for which there exists a solution. We crop it to m_multi_revs
    
    if multi_revs < m_nmax {
        m_nmax = multi_revs;
    }

    // 2.2 We now allocate the memory for the output variables
    let mut m_v1: Vec<[f64; 3]> = Vec::new();
    let mut m_v2: Vec<[f64; 3]> = Vec::new();
    let mut m_iters: Vec<u32> = Vec::new();
    let mut m_x: Vec<f64> = Vec::new();
    for _i in 0..(m_nmax * 2 + 1) {
        m_v1.push([0.; 3]);
        m_v2.push([0.; 3]);
        m_iters.push(0);
        m_x.push(0.);
    }

    // 3 - We may now find all solutions in x,y
    // 3.1 0 rev solution
    // 3.1.1 initial guess
    if t >= t00 { 
        m_x[0] = -(t - t00) / (t - t00 + 4.);
    } else if t <= t1 { 
        m_x[0] = t1 * (t1 - t) / (2. / 5. * (1. - lambda2 * lambda3) * t) + 1.;
    } else {   
        m_x[0] = (t / t00).powf(0.69314718055994529 / (t1 / t00).ln()) - 1.;
    }

    let mut num_iter_tot = 0;
    // 3.1.2 Householder iterations
    let mut m_t = householder_izzo(m_lambda, t, m_x[0], 0, k_p * 1e-3, 50);
    m_x[0] = m_t.0;
    m_iters[0] = m_t.1;

    // 3.2 multi rev solutions
    let mut tmp;
    for i in 0..m_nmax as usize { 
        // 3.2.1 left Householder iterations
        tmp = ((i as f64 * PI + PI) / (8. * t)).powf(2. / 3.);
        m_x[2 * i + 1] = (tmp - 1.) / (tmp + 1.);
        m_t = householder_izzo(m_lambda, t, m_x[2 * i + 1], i as u32, k_p * 1e-3, 50);
        m_x[2 * i + 1] = m_t.0;
        m_iters[2 * i + 1] = m_t.1;

        // 3.2.1 right Householder iterations
        tmp = ((8. * t) / (i as f64 * PI)).powf(2. / 3.);
        m_x[2 * i + 2] = (tmp - 1.) / (tmp + 1.);
        m_t = householder_izzo(m_lambda, t, m_x[2 * i + 2], i as u32, k_p * 1e-3, 50);
        m_x[2 * i + 2] = m_t.0;
        m_iters[2 * i + 2] = m_t.1;
    }

    for j in 0..(2 * m_nmax + 1) as usize {
        num_iter_tot += m_iters[j];
    }

    // 4 - For each found x value we reconstruct the terminal velocities
    let gamma = (MU * m_s / 2.).sqrt();
    let rho = (r1_mag - r2_mag) / m_c;
    let sigma = (1. - rho * rho).sqrt();
    let mut vt;
    let mut vr1;
    let mut vt1;
    let mut vr2;
    let mut vt2;
    let mut y;
    for i in 0..(m_nmax * 2 + 1) as usize {
        y = (1. - lambda2 + lambda2 * m_x[i] * m_x[i]).sqrt();
        vr1 = gamma * ((m_lambda * y - m_x[i]) - rho * (m_lambda * y + m_x[i])) / r1_mag;
        vr2 = -gamma * ((m_lambda * y - m_x[i]) + rho * (m_lambda * y + m_x[i])) / r2_mag;
        vt = gamma * sigma * (y + m_lambda * m_x[i]);
        vt1 = vt / r1_mag;
        vt2 = vt / r2_mag;
        for j in 0..3 as usize {
            m_v1[i][j] = vr1 * ir1[j] + vt1 * it1[j];
        }
        for j in 0..3 as usize {
            m_v2[i][j] = vr2 * ir2[j] + vt2 * it2[j];
        }
    }

    (m_v1, m_v2, num_iter_tot)
}

fn dtdx_izzo(m_lambda: f64, x: f64, t: f64) -> (f64, f64, f64) {
    // solve the value of the derivatives DT, DDT, DDDT 
    // input
    // x:  the value of x
    // T:  non dimensional time-of-flight 

    // output
    // DT:  first derivative
    // DTT:  second Derivative
    // DTTT:  third derivative

    let l2 = m_lambda * m_lambda;
    let l3 = l2 * m_lambda;
    let umx2 = 1. - x * x;
    let y = (1. - l2 * umx2).sqrt();
    let y2 = y * y;
    let y3 = y2 * y;
    let dt   = 1. / umx2 * (3. * t * x - 2. + 2. * l3 * x / y);
    let ddt  = 1. / umx2 * (3. * t + 5. * x * dt + 2. * (1. - l2) * l3 / y3);
    let dddt = 1. / umx2 * (7. * x * ddt + 8. * dt - 6. * (1. - l2) * l2 * l3 * x / y3 / y2);

    (dt, ddt, dddt)
}

fn x2tof_izzo(m_lambda: f64, x: f64, n: u32) -> f64 {
    // solve tof value 
    // input
    //  x:  the value of x
    //  N:  number of multirevs to be computed  
    
    // output
    //  tof:  time of flight
    
    let battin = 0.01;
    let lagrange = 0.2;
    let dist = (x - 1.).abs();
    
    let k = m_lambda * m_lambda;
    let e = x * x - 1.;
    let rho = e.abs();
    let z = (1. + k * e).sqrt();
    
    let tof;
    if dist < lagrange && dist > battin {
        tof = x2tof2_izzo(m_lambda, x, n);
    } else if dist < battin {
        let eta = z - m_lambda * x;
        let s1 = 0.5 * (1. - m_lambda - x * eta);
        let q = hypergeometric_izzo(s1, 1e-11) * 4.0 / 3.0;
        tof = (eta * eta * eta * q + 4. * m_lambda * eta) / 2. + n as f64 * PI / (rho * rho.sqrt());
    } else { 
        let y = rho.sqrt();
        let g = x * z - m_lambda * e;
        let d;
        if e < 0. {
            let l = g.acos();
            d = n as f64 * PI + l;
        } else { 
            let f = y * (z - m_lambda * x);
            d = (f + g).ln();
        }
        tof = (x - m_lambda * z - d / y) / e;
    }
    tof
}

fn x2tof2_izzo(m_lambda: f64, x: f64, n: u32) -> f64 {
    // solve tof value by using Lagrange tof expression
    // input
    // x:  the value of x
    // N:   number of multirevs to be computed   

    // output
    // tof:  time of flight

    let a = 1.0 / (1.0 - x * x);
    let alpha;
    let mut beta;
    let tof;

    if a > 0. {
        // ellipse    
        alpha = 2. * x.acos();
        beta = 2. * (m_lambda * m_lambda / a).sqrt().asin();
        if m_lambda < 0. {
            beta = -beta;
        }
        tof = (a * a.sqrt() * ((alpha - alpha.sin()) - (beta - beta.sin()) + 2. * PI * n as f64)) / 2.;
    } else { 
         alpha = 2. * x.acosh();
         beta = 2. * (-m_lambda * m_lambda / a).sqrt().asinh();
        if m_lambda < 0. { 
            beta = -beta;
        }
        tof = -a * (-a).sqrt() * ((beta - beta.sinh()) - (alpha - alpha.sinh())) / 2.;   
    }      
    tof
}

fn hypergeometric_izzo(z: f64, tol :f64) -> f64 {
    // Gaussian or ordinary hypergeometric function

    let mut cj = 1.;
    let mut sj = 1.;
    let mut err = 1.;
    let mut cj1;
    let mut sj1;
    let mut j = 0.;
    while err > tol {
        cj1 = cj * (3. + j) * (1. + j) / (2.5 + j) * z / (j + 1.);
        sj1 = sj + cj1;
        err = cj1.abs();
        cj = cj1;
        sj = sj1;
        j += 1.;
    }
    sj
}

fn householder_izzo(m_lambda: f64, t: f64, x0: f64, n: u32, eps: f64, iter_max: u32) -> (f64, u32) {
    // solve x value by using Householder iteration
    // input
    // T:  non dimensional time-of-flight
    // x0:  initial guess 
    // N:  number of multirevs to be computed  
    // eps:  Error precision for TolF
    // iter_max:  Maximum number of iterations

    // output
    // x:  solution for x
    // it:  Number of actual iterations

    let mut it = 0;
    let mut err = 1.;
    let mut x0_ = x0;
    let e3 = 1e-6;
    let mut xnew;
    let mut tof;
    let mut dt;
    let mut ddt;
    let mut dddt;
    let mut dt2;
    let mut dt3;
    let mut delta: f64 = 10.;
    while (err > e3 || delta.abs() > eps) && (it < iter_max) {
        tof = x2tof_izzo(m_lambda, x0_, n);
        dt3 = dtdx_izzo(m_lambda, x0_, tof);
        dt = dt3.0;
        ddt = dt3.1;
        dddt = dt3.2;
        delta = tof - t;
        dt2 = dt * dt;
        xnew = x0_ - delta * (dt2 - delta * ddt / 2.) / (dt * (dt2 - delta * ddt) + dddt * delta * delta / 6.);
        err = (x0_ - xnew).abs();
        x0_ = xnew;
        it += 1;
    }

    (x0_, it)
}