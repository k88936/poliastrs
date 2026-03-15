use nalgebra::{Matrix3, Vector3};

use crate::twobody::states::CartesianState;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassicalElements {
    pub p_km: f64,
    pub ecc: f64,
    pub inc_rad: f64,
    pub raan_rad: f64,
    pub argp_rad: f64,
    pub nu_rad: f64,
}

pub fn coe2rv(mu_km3_s2: f64, coe: ClassicalElements) -> CartesianState {
    let p = coe.p_km;
    let e = coe.ecc;
    let nu = coe.nu_rad;

    let r_pqw = Vector3::new(
        p * nu.cos() / (1.0 + e * nu.cos()),
        p * nu.sin() / (1.0 + e * nu.cos()),
        0.0,
    );
    let v_pqw = Vector3::new(
        -(mu_km3_s2 / p).sqrt() * nu.sin(),
        (mu_km3_s2 / p).sqrt() * (e + nu.cos()),
        0.0,
    );

    let (raan, inc, argp) = (coe.raan_rad, coe.inc_rad, coe.argp_rad);
    let rot = rotation_313(raan, inc, argp);
    CartesianState::new(rot * r_pqw, rot * v_pqw)
}

pub fn rv2coe(mu_km3_s2: f64, state: &CartesianState) -> ClassicalElements {
    let r = state.r_km;
    let v = state.v_km_s;
    let r_norm = r.norm();
    let v_norm = v.norm();
    let h = r.cross(&v);
    let h_norm = h.norm();
    let k_hat = Vector3::new(0.0, 0.0, 1.0);
    let n = k_hat.cross(&h);
    let n_norm = n.norm();
    let e_vec = ((v_norm * v_norm - mu_km3_s2 / r_norm) * r - r.dot(&v) * v) / mu_km3_s2;
    let ecc = e_vec.norm();
    let p_km = h_norm * h_norm / mu_km3_s2;
    let inc_rad = (h.z / h_norm).acos();

    let raan_rad = if n_norm > 1e-14 {
        n.y.atan2(n.x).rem_euclid(std::f64::consts::TAU)
    } else {
        0.0
    };

    let argp_rad = if n_norm > 1e-14 && ecc > 1e-14 {
        let c = (n.dot(&e_vec) / (n_norm * ecc)).clamp(-1.0, 1.0);
        let mut w = c.acos();
        if e_vec.z < 0.0 {
            w = std::f64::consts::TAU - w;
        }
        w
    } else {
        0.0
    };

    let nu_rad = if ecc > 1e-14 {
        let c = (e_vec.dot(&r) / (ecc * r_norm)).clamp(-1.0, 1.0);
        let mut nu = c.acos();
        if r.dot(&v) < 0.0 {
            nu = std::f64::consts::TAU - nu;
        }
        nu
    } else {
        let c = if n_norm > 1e-14 {
            (n.dot(&r) / (n_norm * r_norm)).clamp(-1.0, 1.0)
        } else {
            (r.x / r_norm).clamp(-1.0, 1.0)
        };
        let mut nu = c.acos();
        if r.z < 0.0 {
            nu = std::f64::consts::TAU - nu;
        }
        nu
    };

    ClassicalElements {
        p_km,
        ecc,
        inc_rad,
        raan_rad,
        argp_rad,
        nu_rad,
    }
}

fn rotation_313(raan: f64, inc: f64, argp: f64) -> Matrix3<f64> {
    let (co, so) = (raan.cos(), raan.sin());
    let (ci, si) = (inc.cos(), inc.sin());
    let (cw, sw) = (argp.cos(), argp.sin());
    Matrix3::new(
        co * cw - so * sw * ci,
        -co * sw - so * cw * ci,
        so * si,
        so * cw + co * sw * ci,
        -so * sw + co * cw * ci,
        -co * si,
        sw * si,
        cw * si,
        ci,
    )
}
