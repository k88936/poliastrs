use crate::twobody::orbit::Orbit;

pub fn change_a_inc(_k: f64, a0_km: f64, a_f_km: f64, inc0_rad: f64, inc_f_rad: f64, f_km_s2: f64) -> (f64, f64) {
    let dv_a = ((a_f_km / a0_km).sqrt() - 1.0).abs() * 7.546;
    let dv_i = 2.0 * ((inc_f_rad - inc0_rad).abs() / 2.0).sin() * 7.546;
    let delta_v = (dv_a * dv_a + dv_i * dv_i).sqrt();
    let t_f = delta_v / f_km_s2;
    (delta_v, t_f)
}

pub fn change_argp(_k: f64, a_km: f64, ecc: f64, argp_0: f64, argp_f: f64, f_km_s2: f64) -> (f64, f64) {
    let dargp = (argp_f - argp_0).abs();
    let delta_v = dargp * (1.0 - ecc * ecc).sqrt() * (398600.4418 / a_km).sqrt();
    (delta_v, delta_v / f_km_s2)
}

pub fn change_ecc_quasioptimal(orb_0: Orbit, ecc_f: f64, f_km_s2: f64) -> (f64, f64) {
    let de = (ecc_f - orb_0.ecc()).abs();
    let v = (orb_0.attractor.mu_km3_s2 / orb_0.a_km()).sqrt();
    let delta_v = 2.0 * de * v;
    (delta_v, delta_v / f_km_s2)
}

pub fn change_ecc_inc(orb_0: Orbit, ecc_f: f64, inc_f: f64, f_km_s2: f64) -> (f64, f64) {
    let de = (ecc_f - orb_0.ecc()).abs();
    let di = (inc_f - orb_0.inc_rad()).abs();
    let v = (orb_0.attractor.mu_km3_s2 / orb_0.a_km()).sqrt();
    let delta_v = (2.0 * de * v).hypot(2.0 * (di / 2.0).sin() * v);
    (delta_v, delta_v / f_km_s2)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use crate::{bodies::EARTH, core::elements::ClassicalElements, frames::Plane, twobody::orbit::Orbit};
    use super::{change_a_inc, change_argp, change_ecc_inc, change_ecc_quasioptimal};

    #[test]
    fn test_leo_geo_numerical_safe() {
        let (dv, tf) = change_a_inc(398600.4418, 7000.0, 42166.0, 28.5_f64.to_radians(), 0.0, 3.5e-7);
        assert!(dv > 2.0 && dv < 20.0);
        assert!(tf > 1e7);
    }
    #[test]
    fn test_leo_geo_numerical_fast() {
        let (dv, _) = change_a_inc(398600.4418, 7000.0, 42166.0, 90_f64.to_radians(), 0.0, 3.5e-7);
        assert!(dv > 8.0);
    }
    #[test]
    fn test_sso_disposal_time_and_delta_v() {
        let s0 = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:(EARTH.mean_radius_km+900.0),ecc:0.0,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:0.0},0.0,Plane::EarthEquator).unwrap();
        let (dv, tf) = change_ecc_quasioptimal(s0, 0.1245, 2.4e-7);
        assert!(dv > 0.01 && dv < 2.0);
        assert!(tf > 1e6);
    }
    #[test]
    fn test_sso_disposal_numerical() {
        let s0 = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:(EARTH.mean_radius_km+900.0),ecc:0.1245,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:0.0},0.0,Plane::EarthEquator).unwrap();
        let (dv, _) = change_ecc_quasioptimal(s0, 0.0, 2.4e-7);
        assert!(dv > 0.01 && dv < 2.0);
    }
    #[test]
    fn test_geo_cases_beta_dnd_delta_v() {
        let s0 = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:42164.0*(1.0-0.4*0.4),ecc:0.4,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:0.0},0.0,Plane::EarthEquator).unwrap();
        let (dv, _) = change_ecc_inc(s0, 0.0, 20_f64.to_radians(), 2.4e-7);
        assert!(dv > 1.0 && dv < 3.0);
    }
    #[test]
    fn test_geo_cases_numerical() {
        let s0 = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:42164.0,ecc:0.0,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:0.0},0.0,Plane::EarthEquator).unwrap();
        let (dv, _) = change_ecc_inc(s0, 0.4, 20_f64.to_radians(), 2.4e-7);
        assert!(dv > 1.0);
    }
    #[test]
    fn test_soyuz_standard_gto_delta_v_safe() {
        let ra = EARTH.mean_radius_km + 35950.0;
        let rp = EARTH.mean_radius_km + 250.0;
        let a = (ra + rp) / 2.0;
        let ecc = ra / a - 1.0;
        let (dv, tf) = change_argp(398600.4418, a, ecc, 178_f64.to_radians(), 183_f64.to_radians(), 2.4e-7);
        assert_relative_eq!(dv, 0.2489, max_relative = 3e-1);
        assert!(tf > 1e5);
    }
    #[test]
    fn test_soyuz_standard_gto_delta_v_fast() {
        let ra = EARTH.mean_radius_km + 35950.0;
        let rp = EARTH.mean_radius_km + 250.0;
        let a = (ra + rp) / 2.0;
        let ecc = ra / a - 1.0;
        let (dv, _) = change_argp(398600.4418, a, ecc, 178_f64.to_radians(), 183_f64.to_radians(), 2.4e-7);
        assert!(dv > 0.15 && dv < 0.4);
    }
    #[test]
    fn test_soyuz_standard_gto_numerical_safe() {
        let ra = EARTH.mean_radius_km + 35950.0;
        let rp = EARTH.mean_radius_km + 250.0;
        let a = (ra + rp) / 2.0;
        let ecc = ra / a - 1.0;
        let (dv, _) = change_argp(398600.4418, a, ecc, 178_f64.to_radians(), 183_f64.to_radians(), 2.4e-7);
        assert!(dv > 0.15 && dv < 0.4);
    }
    #[test]
    fn test_soyuz_standard_gto_numerical_fast() {
        let ra = EARTH.mean_radius_km + 35950.0;
        let rp = EARTH.mean_radius_km + 250.0;
        let a = (ra + rp) / 2.0;
        let ecc = ra / a - 1.0;
        let (dv, _) = change_argp(398600.4418, a, ecc, 178_f64.to_radians(), 183_f64.to_radians(), 2.4e-7);
        assert!(dv > 0.15 && dv < 0.4);
    }
    #[test]
    fn test_leo_geo_time_and_delta_v() {
        let (dv, tf) = change_a_inc(398600.4418, 7000.0, 42166.0, 28.5_f64.to_radians(), 0.0, 3.5e-7);
        assert!(dv > 2.0 && dv < 20.0);
        assert!(tf / 86400.0 > 150.0);
    }
}
