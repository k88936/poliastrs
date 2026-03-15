use crate::core::angles::{e_to_m, f_to_m, m_to_e, m_to_f, nu_to_e, nu_to_f};
use crate::core::elements::{ClassicalElements, coe2rv, rv2coe};
use crate::twobody::states::CartesianState;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ValladoError {
    ParabolicUnsupported,
}

pub fn propagate_vallado(
    mu_km3_s2: f64,
    state: &CartesianState,
    dt_seconds: f64,
) -> Result<CartesianState, ValladoError> {
    let coe = rv2coe(mu_km3_s2, state);
    if (coe.ecc - 1.0).abs() < 1e-10 {
        return Err(ValladoError::ParabolicUnsupported);
    }
    let (new_nu, p_km) = if coe.ecc < 1.0 {
        let a = coe.p_km / (1.0 - coe.ecc * coe.ecc);
        let n = (mu_km3_s2 / (a * a * a)).sqrt();
        let m0 = e_to_m(nu_to_e(coe.nu_rad, coe.ecc), coe.ecc);
        let m = m0 + n * dt_seconds;
        let e = m_to_e(m, coe.ecc);
        let nu = crate::core::angles::e_to_nu(e, coe.ecc);
        (nu, coe.p_km)
    } else {
        let a_abs = coe.p_km / (coe.ecc * coe.ecc - 1.0);
        let n = (mu_km3_s2 / (a_abs * a_abs * a_abs)).sqrt();
        let m0 = f_to_m(nu_to_f(coe.nu_rad, coe.ecc), coe.ecc);
        let m = m0 + n * dt_seconds;
        let f = m_to_f(m, coe.ecc);
        let nu = crate::core::angles::f_to_nu(f, coe.ecc);
        (nu, coe.p_km)
    };
    Ok(coe2rv(
        mu_km3_s2,
        ClassicalElements {
            nu_rad: new_nu,
            p_km,
            ecc: coe.ecc,
            inc_rad: coe.inc_rad,
            raan_rad: coe.raan_rad,
            argp_rad: coe.argp_rad,
        },
    ))
}
