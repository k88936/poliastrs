#[cfg(test)]
mod tests {
    use crate::examples::iss;
    use crate::twobody::propagation::farnocchia::propagate_two_body;
    use crate::twobody::orbit::Orbit;
    use approx::assert_relative_eq;

    #[test]
    fn test_propagate_farnocchia_iss_period() {
        let orbit = iss();
        let period = orbit.period_seconds().unwrap();
        let mu = orbit.attractor.mu_km3_s2;
        
        let state = orbit.state;
        let final_state_res = propagate_two_body(mu, &state, period);
        assert!(final_state_res.is_ok());
        let final_state = final_state_res.unwrap();
        
        // Convert back to orbit to check nu
        let final_orbit = Orbit::from_vectors_at(
            orbit.attractor,
            [final_state.r_km.x, final_state.r_km.y, final_state.r_km.z],
            [final_state.v_km_s.x, final_state.v_km_s.y, final_state.v_km_s.z],
            orbit.epoch_tdb_seconds + period, // epoch updated
            orbit.plane,
        );
        
        // nu should be same as initial nu (modulo 2pi)
        let nu_initial = orbit.classical().nu_rad.to_degrees();
        let nu_final = final_orbit.classical().nu_rad.to_degrees();
        
        // Normalize angles to [0, 360) or handles wraparound
        let diff = (nu_final - nu_initial).abs() % 360.0;
        let diff = if diff > 180.0 { 360.0 - diff } else { diff };
        
        assert_relative_eq!(diff, 0.0, epsilon = 1e-4);
        
        // Also check Cartesian state matches
        assert_relative_eq!(final_state.r_km, state.r_km, epsilon = 1e-6);
        assert_relative_eq!(final_state.v_km_s, state.v_km_s, epsilon = 1e-6);
    }
}
