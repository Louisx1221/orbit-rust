// mod math;
use crate::math::*;

pub const MU: f32 = 3.986004118e14;

pub fn solve_kepler(ma: f32, ecc: f32) -> f32 {
    // 初始参数
    // 设置容错度
    let tol = 1e-8;
    // 设置最大迭代次数
    let iter_max = 100;

    // 设置偏近点角E初值
    let mut ea: f32 = fmod(ma, TWO_PI);
    if ea < PI {
        ea += ecc / 2.0;
    } else {
        ea -= ecc / 2.0;
    }

    // 牛顿迭代
    // 迭代次数
    let mut iter = 1;
    // 斜率
    let mut ratio: f32 = 1.0;
    while (ratio.abs() > tol) || (iter < iter_max) {
        ratio = (ea - ecc * ea.sin() - ma) / (1.0 - ecc * ea.cos());
        ea -= ratio;
        iter += 1;
    }

    ea
}

# [derive(Debug)]
pub struct Orb {
    pub r: [f32; 3],
    pub v: [f32; 3],
    pub coe: [f32; 6],
    pub jd: f32,
}

impl Orb {
    pub fn rv2coe(&mut self) {
        // 位置及速度大小
        let r_mag = norm(self.r);
        let v_mag = norm(self.v);

        // 位置及速度单位向量
        let mut r_hat: [f32; 3] = self.r;
        for r_hat_i in r_hat.iter_mut() {
            *r_hat_i /= r_mag;
        }
 
        // 角动量
        let hv = cross(self.r, self.v);
        let hv_mag = norm(hv);
        let mut h_hat: [f32; 3] = hv;
        for h_hat_i in h_hat.iter_mut() {
            *h_hat_i /= hv_mag;
        }

        // 偏心率矢量
        let mut v_tmp: [f32; 3] = self.v;
        for v_tmp_i in v_tmp.iter_mut() {
            *v_tmp_i /= MU;
        }
        let mut ecc_vec = cross(v_tmp, hv);
        for (i, ecc_vec_i) in ecc_vec.iter_mut().enumerate() {
            *ecc_vec_i -= r_hat[i];
        }

        // 半长轴
        let sma = 1.0 / (2.0 / r_mag - v_mag * v_mag / MU);

        let p = if h_hat[0] == 0.0 {
            0.0
        } else {
            h_hat[0] / (1.0 + h_hat[2])
        };
        let q = if h_hat[1] == 0.0 {
            0.0
        } else {
            -h_hat[1] / (1.0 + h_hat[2])
        };

        let const1 = 1.0 / (1.0 + p * p + q * q);
        let mut f_hat: [f32; 3] = [0.0; 3];
        let mut g_hat: [f32; 3] = [0.0; 3];
        f_hat[0] = const1 * (1.0 - p * p + q * q);
        f_hat[1] = const1 * 2.0 * p * q;
        f_hat[2] = -const1 * 2.0 * p;
        g_hat[0] = const1 * 2.0 * p * q;
        g_hat[1] = const1 * (1.0 + p * p - q * q);
        g_hat[2] = const1 * 2.0 * q;
        let h = dot(ecc_vec, g_hat);
        let xk = dot(ecc_vec, f_hat);
        let x1 = dot(self.r, f_hat);
        let y1 = dot(self.r, g_hat);

        // 偏心率
        let ecc = (h * h + xk * xk).sqrt();

        // 轨道倾角
        let mut inc = fmod(2.0 * (p * p + q * q).sqrt().atan(), TWO_PI); 

        // 真实纬度
        let lambda_t = fmod(atan2(y1, x1), TWO_PI);

        // 升交点赤经
        let mut raan = if inc > 1e-8 {
            fmod(atan2(p, q), TWO_PI)
        } else {
            0.0
        };

        // 近地点幅角
        let mut aop = if ecc > 1e-8 {
            fmod(fmod(atan2(h, xk), TWO_PI) - raan, TWO_PI)
        } else {
            0.0
        };

        // 真近点角
        let mut ta = fmod(lambda_t - raan - aop, TWO_PI); 

        // 奇异值
        if h_hat[2] == -1.0 {
            inc  = PI;
            raan = PI;
            aop  = fmod(-fmod(atan2(h, xk), TWO_PI) + raan, TWO_PI);
            ta   = fmod(lambda_t - raan - aop, TWO_PI);
        }

        // 经典轨道要素
        self.coe[0] = sma;
        self.coe[1] = ecc;
        self.coe[2] = inc;
        self.coe[3] = raan;
        self.coe[4] = aop;
        self.coe[5] = ta;
    }
}