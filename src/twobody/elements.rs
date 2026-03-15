pub fn circular_velocity(mu_km3_s2: f64, a_km: f64) -> f64 {
    (mu_km3_s2 / a_km).sqrt()
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::circular_velocity;

    #[test]
    fn simple_circular_velocity() {
        let v = circular_velocity(398600.0, 7000.0);
        assert_relative_eq!(v, 7.5460491, max_relative = 1e-7);
    }
}
