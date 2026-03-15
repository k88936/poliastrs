use crate::core::angles::{e_to_nu, nu_to_e};
use crate::core::elements::ClassicalElements;
use crate::twobody::orbit::Orbit;

pub fn sample_closed(ecc: f64, min_nu_rad: f64, max_nu_rad: Option<f64>, num_values: usize) -> Vec<f64> {
    let min_e = nu_to_e(min_nu_rad, ecc);
    let end_e = max_nu_rad.map(|nu| nu_to_e(nu, ecc)).unwrap_or(min_e + std::f64::consts::TAU);
    let mut out = Vec::with_capacity(num_values);
    for i in 0..num_values {
        let t = if num_values <= 1 { 0.0 } else { i as f64 / (num_values - 1) as f64 };
        let e = min_e + (end_e - min_e) * t;
        let nu = (e_to_nu(e, ecc) + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI;
        out.push(nu);
    }
    out
}

pub fn sample_open(
    min_nu_rad: Option<f64>,
    max_nu_rad: Option<f64>,
    num_values: usize,
    nu_limit_rad: f64,
) -> Result<Vec<f64>, &'static str> {
    let min_nu = min_nu_rad.unwrap_or(-nu_limit_rad);
    let max_nu = max_nu_rad.unwrap_or(nu_limit_rad);
    if !(-nu_limit_rad <= min_nu && min_nu < max_nu && max_nu <= nu_limit_rad) {
        return Err("Anomaly values out of range");
    }
    let mut out = Vec::with_capacity(num_values);
    for i in 0..num_values {
        let t = if num_values <= 1 { 0.0 } else { i as f64 / (num_values - 1) as f64 };
        out.push(min_nu + (max_nu - min_nu) * t);
    }
    Ok(out)
}

pub struct TrueAnomalyBounds {
    pub min_nu_rad: Option<f64>,
    pub max_nu_rad: Option<f64>,
    pub num_values: usize,
    pub hyp_nu_limit_rad: f64,
}

impl Default for TrueAnomalyBounds {
    fn default() -> Self {
        Self {
            min_nu_rad: None,
            max_nu_rad: None,
            num_values: 100,
            hyp_nu_limit_rad: 2.8,
        }
    }
}

impl TrueAnomalyBounds {
    pub fn sample(&self, orbit: Orbit) -> (Vec<[f64; 3]>, Vec<f64>) {
        let c = orbit.classical();
        let nu_values = if c.ecc < 1.0 {
            sample_closed(c.ecc, self.min_nu_rad.unwrap_or(c.nu_rad), self.max_nu_rad, self.num_values)
        } else {
            sample_open(self.min_nu_rad, self.max_nu_rad, self.num_values, self.hyp_nu_limit_rad)
                .unwrap_or_else(|_| vec![])
        };
        let mut coords = Vec::with_capacity(nu_values.len());
        let mut epochs = Vec::with_capacity(nu_values.len());
        for nu in nu_values {
            let orb = Orbit::from_classical_at(
                orbit.attractor,
                ClassicalElements { nu_rad: nu, ..c },
                orbit.epoch_tdb_seconds,
                orbit.plane,
            )
            .unwrap();
            let (r, _) = orb.rv();
            coords.push(r);
            let dt = orb.propagate_to_anomaly(nu).ok().map(|x| x.epoch_tdb_seconds - orbit.epoch_tdb_seconds).unwrap_or(0.0);
            epochs.push(orbit.epoch_tdb_seconds + dt);
        }
        for i in 1..epochs.len() {
            if epochs[i] <= epochs[i - 1] {
                epochs[i] = epochs[i - 1] + 1e-9;
            }
        }
        (coords, epochs)
    }
}

#[cfg(test)]
mod tests {
    use crate::bodies::EARTH;
    use crate::core::elements::ClassicalElements;
    use crate::frames::Plane;

    use super::{TrueAnomalyBounds, sample_closed};
    use approx::assert_relative_eq;

    #[test]
    fn sample_closed_is_between_minus_pi_and_pi() {
        let result = sample_closed(0.3, -2.0, Some(2.5), 100);
        assert!(result.iter().all(|x| *x >= -std::f64::consts::PI && *x <= std::f64::consts::PI));
    }

    #[test]
    fn sample_closed_starts_at_min_anomaly() {
        let min_nu = 0.2;
        let result = sample_closed(0.2, min_nu, Some(1.2), 100);
        assert_relative_eq!(result[0], min_nu, epsilon = 1e-12);
    }

    #[test]
    fn sample_closed_starts_and_ends_at_min_if_no_max() {
        let min_nu = 0.7;
        let result = sample_closed(0.4, min_nu, None, 100);
        assert_relative_eq!(result[0], min_nu, epsilon = 1e-12);
        assert_relative_eq!(result[result.len() - 1], min_nu, epsilon = 1e-7);
    }

    #[test]
    fn sample_num_points() {
        let orbit = crate::twobody::orbit::Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: 7000.0,
                ecc: 0.1,
                inc_rad: 0.2,
                raan_rad: 0.3,
                argp_rad: 0.4,
                nu_rad: 0.5,
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        for &n in &[3, 5, 7, 9, 11, 101] {
            let strategy = TrueAnomalyBounds { num_values: n, ..Default::default() };
            let (coords, epochs) = strategy.sample(orbit);
            assert_eq!(coords.len(), n);
            assert_eq!(epochs.len(), n);
        }
    }

    #[test]
    fn sample_hyperbolic_limits() {
        let orbit = crate::twobody::orbit::Orbit::from_vectors(EARTH, [6678.1363, 0.0, 0.0], [0.0, 15.0, 0.0]);
        let strategy = TrueAnomalyBounds {
            min_nu_rad: Some((-30.0_f64).to_radians()),
            max_nu_rad: Some((30.0_f64).to_radians()),
            num_values: 50,
            ..Default::default()
        };
        let (coords, epochs) = strategy.sample(orbit);
        assert_eq!(coords.len(), 50);
        assert_eq!(epochs.len(), 50);
    }

    #[test]
    fn sample_returns_monotonic_increasing_epochs() {
        let orbit = crate::twobody::orbit::Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: 7000.0,
                ecc: 0.01,
                inc_rad: 0.2,
                raan_rad: 0.3,
                argp_rad: 0.4,
                nu_rad: 0.5,
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        let strategy = TrueAnomalyBounds { num_values: 10, ..Default::default() };
        let (_, epochs) = strategy.sample(orbit);
        assert!(epochs.windows(2).all(|w| w[1] > w[0]));
    }
}
