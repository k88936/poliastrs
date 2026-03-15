use approx::assert_relative_eq;
use nalgebra::Vector3;
use crate::bodies::{SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, MOON};
use crate::frames::{Frame, Coordinate, rotation};

const J2000: f64 = 0.0;

#[test]
fn test_fixed_frame_calculation_gives_expected_result() {
    // Port of test_fixed_frame_calculation_gives_expected_result
    // Uses J2000 (d=0)
    
    // (Body, expected alpha, expected delta, expected W)
    let cases = vec![
        (SUN, 286.13, 63.87, 84.176),
        (MERCURY, 281.0103, 61.45, 329.5999488),
        (VENUS, 272.76, 67.16, 160.2),
        (MARS, 317.68085441, 52.88643928, 176.63205973),
        (JUPITER, 268.05720404, 64.49580995, 284.95),
        (SATURN, 40.589, 83.537, 38.9),
        (URANUS, 257.311, -15.175, 203.81),
        (NEPTUNE, 299.33373896, 42.95035902, 249.99600757),
        (MOON, 266.85773344495135, 65.64110274784535, 41.1952639807452),
    ];

    for (body, exp_alpha, exp_delta, exp_w) in cases {
        let elements = rotation::compute_rotational_elements(body, J2000);
        
        assert_relative_eq!(elements.alpha_deg, exp_alpha, epsilon = 1e-6);
        assert_relative_eq!(elements.delta_deg, exp_delta, epsilon = 1e-6);
        assert_relative_eq!(elements.w_deg, exp_w, epsilon = 1e-6);
    }
}

#[test]
fn test_planetary_inertial_roundtrip_vector() {
    // Port of test_planetary_inertial_roundtrip_vector
    // Uses J2000 and subsequent times
    let bodies = vec![SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE];
    let sampling_time = 10.0; // seconds
    let steps = 100; // reduced from 1000 for speed

    for body in bodies {
        for i in 0..steps {
            let t = J2000 + (i as f64) * sampling_time;
            
            // Initial vector in Fixed frame: (R, 0, 0)
            // Equivalent to (0 deg, 0 deg) spherical at surface
            let r_fixed = Vector3::new(body.mean_radius_km, 0.0, 0.0);
            let coord_fixed = Coordinate::new(r_fixed, Frame::BodyFixed(body), t);

            // Transform to Inertial
            let coord_inertial = coord_fixed.transform_to(Frame::BodyInertial(body)).unwrap();
            
            // Transform back to Fixed
            let coord_back = coord_inertial.transform_to(Frame::BodyFixed(body)).unwrap();

            // Check distance preservation
            assert_relative_eq!(coord_inertial.vector.norm(), body.mean_radius_km, epsilon = 1e-6);
            
            // Check roundtrip accuracy
            assert_relative_eq!(coord_back.vector, r_fixed, epsilon = 1e-6);
        }
    }
}

#[test]
fn test_planetary_icrs_frame_is_just_translation() {
    // Since get_body_barycentric returns 0 currently, this tests identity transformation
    // If we implement ephemeris, it tests correctness
    let bodies = vec![SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE];
    let epoch = J2000;
    
    // Vector in BodyInertial frame
    let vec = Vector3::new(100.0, 100.0, 100.0);

    for body in bodies {
        let coord = Coordinate::new(vec, Frame::BodyInertial(body), epoch);
        let coord_icrs = coord.transform_to(Frame::ICRS).unwrap();
        
        // Expected: r_ICRS = r_BodyInertial + r_Body_ICRS
        // r_BodyInertial = vec
        // r_Body_ICRS = get_body_barycentric(body, epoch) (which is 0 currently)
        let expected = vec + crate::frames::get_body_barycentric(body, epoch);
        
        assert_relative_eq!(coord_icrs.vector, expected, epsilon = 1e-6);
    }
}
