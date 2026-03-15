use std::f64::consts::PI;

pub fn e_to_m(e_anomaly: f64, ecc: f64) -> f64 {
    e_anomaly - ecc * e_anomaly.sin()
}

pub fn m_to_e(mean_anomaly: f64, ecc: f64) -> f64 {
    let mut e = mean_anomaly;
    for _ in 0..100 {
        let f = e - ecc * e.sin() - mean_anomaly;
        let fp = 1.0 - ecc * e.cos();
        let de = f / fp;
        e -= de;
        if de.abs() < 1e-13 {
            break;
        }
    }
    e
}

pub fn nu_to_e(nu: f64, ecc: f64) -> f64 {
    let beta = ((1.0 - ecc) / (1.0 + ecc)).sqrt();
    2.0 * (beta * (nu / 2.0).tan()).atan()
}

pub fn e_to_nu(e_anomaly: f64, ecc: f64) -> f64 {
    let beta = ((1.0 + ecc) / (1.0 - ecc)).sqrt();
    2.0 * (beta * (e_anomaly / 2.0).tan()).atan()
}

pub fn f_to_m(f_anomaly: f64, ecc: f64) -> f64 {
    ecc * f_anomaly.sinh() - f_anomaly
}

pub fn m_to_f(mean_anomaly: f64, ecc: f64) -> f64 {
    let mut f = (mean_anomaly / ecc.max(1.01)).asinh();
    for _ in 0..100 {
        let h = ecc * f.sinh() - f - mean_anomaly;
        let hp = ecc * f.cosh() - 1.0;
        let df = h / hp;
        f -= df;
        if df.abs() < 1e-13 {
            break;
        }
    }
    f
}

pub fn nu_to_f(nu: f64, ecc: f64) -> f64 {
    let s = ((ecc - 1.0) / (ecc + 1.0)).sqrt() * (nu / 2.0).tan();
    2.0 * s.atanh()
}

pub fn f_to_nu(f_anomaly: f64, ecc: f64) -> f64 {
    let s = ((ecc + 1.0) / (ecc - 1.0)).sqrt() * (f_anomaly / 2.0).tanh();
    2.0 * s.atan()
}

pub fn wrap_to_pi(angle: f64) -> f64 {
    (angle + PI).rem_euclid(2.0 * PI) - PI
}
