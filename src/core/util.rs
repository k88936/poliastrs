use std::f64::consts::PI;
use nalgebra::{Matrix3, Vector3, Rotation3};

pub fn rotation_matrix(angle_rad: f64, axis: usize) -> Matrix3<f64> {
    let axis_vec = match axis {
        0 => Vector3::x_axis(),
        1 => Vector3::y_axis(),
        2 => Vector3::z_axis(),
        _ => panic!("Invalid axis: {}", axis),
    };
    Rotation3::from_axis_angle(&axis_vec, angle_rad).into_inner()
}

pub fn spherical_to_cartesian(r_colat_lon: Vector3<f64>) -> Vector3<f64> {
    // Input: [r, colat, lon] (radians)
    // colat is colatitude (angle from Z axis, 0 to pi)
    // lon is azimuthal angle (0 to 2pi)
    
    let r = r_colat_lon[0];
    let colat = r_colat_lon[1];
    let lon = r_colat_lon[2];
    
    let x = r * colat.sin() * lon.cos();
    let y = r * colat.sin() * lon.sin();
    let z = r * colat.cos();
    
    Vector3::new(x, y, z)
}

pub fn planetocentric_to_altaz(theta: f64, phi: f64) -> Matrix3<f64> {
    // Transformation matrix from Planetocentric to AltAz
    // theta: Local Sidereal Time
    // phi: Planetodetic latitude
    
    let sin_theta = theta.sin();
    let cos_theta = theta.cos();
    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    
    Matrix3::new(
        -sin_theta, cos_theta, 0.0,
        -sin_phi * cos_theta, -sin_phi * sin_theta, cos_phi,
        cos_phi * cos_theta, cos_phi * sin_theta, sin_phi,
    )
}

pub fn alinspace(start: f64, end: Option<f64>, num: usize) -> Vec<f64> {
    let mut result = Vec::with_capacity(num);
    
    if num == 0 {
        return result;
    }
    
    let start_val = start;
    let end_val = if let Some(e) = end {
        if e <= start {
             // Wrap to next forward angle
             let diff = (e - start).rem_euclid(2.0 * PI);
             if diff.abs() < 1e-12 && e != start {
                 start + 2.0 * PI
             } else {
                 start + diff
             }
        } else {
            e
        }
    } else {
        start + 2.0 * PI
    };
    
    if num == 1 {
        result.push(start_val);
        return result;
    }
    
    let step = (end_val - start_val) / ((num - 1) as f64);
    for i in 0..num {
        result.push(start_val + (i as f64) * step);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_rotation_matrix_axes() {
        let angle = PI / 2.0;
        // Rx(90) = [[1,0,0],[0,0,-1],[0,1,0]]
        let rx = rotation_matrix(angle, 0);
        assert_relative_eq!(rx[(0,0)], 1.0);
        assert_relative_eq!(rx[(1,1)], 0.0);
        assert_relative_eq!(rx[(1,2)], -1.0);
        
        // Rz(90) = [[0,-1,0],[1,0,0],[0,0,1]]
        let rz = rotation_matrix(angle, 2);
        assert_relative_eq!(rz[(0,1)], -1.0);
        assert_relative_eq!(rz[(1,0)], 1.0);
    }

    #[test]
    fn test_spherical_to_cartesian() {
        // [0.5, pi/4, -pi/4] -> [0.25, -0.25, 0.3535]
        let input = Vector3::new(0.5, PI/4.0, -PI/4.0);
        let res = spherical_to_cartesian(input);
        
        assert_relative_eq!(res.x, 0.25);
        assert_relative_eq!(res.y, -0.25);
        assert_relative_eq!(res.z, 0.35355339059, epsilon = 1e-9);
    }

    #[test]
    fn test_alinspace_always_increasing() {
        let res = alinspace(0.0, Some(-0.1), 10);
        for i in 0..res.len()-1 {
            assert!(res[i+1] >= res[i]);
        }
        // Check end value matches wrapped end
        let end_expected = 2.0 * PI - 0.1;
        assert_relative_eq!(*res.last().unwrap(), end_expected, epsilon=1e-6);
    }

    #[test]
    fn test_alinspace_full_circle() {
        let res = alinspace(0.0, None, 5);
        // 0, pi/2, pi, 3pi/2, 2pi
        assert_relative_eq!(res[0], 0.0);
        assert_relative_eq!(*res.last().unwrap(), 2.0*PI);
    }
}
