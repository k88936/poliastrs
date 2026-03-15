use nalgebra::{Vector3, Matrix3, Rotation3};
use crate::bodies::Body;
use std::f64::consts::FRAC_PI_2;

pub mod rotation;
use rotation::{compute_rotational_elements, RotationalElements};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Frame {
    ICRS,
    GCRS,
    HeliocentricInertial,
    BodyInertial(Body),
    BodyFixed(Body),
    GeocentricSolarEcliptic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Plane {
    EarthEquator,
    EarthEcliptic,
}

#[derive(Debug, Clone, Copy)]
pub struct Coordinate {
    pub vector: Vector3<f64>, // km
    pub frame: Frame,
    pub epoch: f64, // TDB seconds from J2000
}

impl Coordinate {
    pub fn new(vector: Vector3<f64>, frame: Frame, epoch: f64) -> Self {
        Self { vector, frame, epoch }
    }

    pub fn transform_to(self, target_frame: Frame) -> Result<Self, String> {
        if self.frame == target_frame {
            return Ok(self);
        }

        match (self.frame, target_frame) {
            // BodyFixed -> BodyInertial (Same Body)
            (Frame::BodyFixed(b1), Frame::BodyInertial(b2)) if b1.name == b2.name => {
                let rot_elements = compute_rotational_elements(b1, self.epoch);
                // M transforms Inertial -> Fixed
                // So Fixed -> Inertial uses M^T
                let rot_matrix = rotation_matrix_inertial_to_fixed(rot_elements);
                let v_inertial = rot_matrix.transpose() * self.vector;
                Ok(Coordinate::new(v_inertial, target_frame, self.epoch))
            },
            // BodyInertial -> BodyFixed (Same Body)
            (Frame::BodyInertial(b1), Frame::BodyFixed(b2)) if b1.name == b2.name => {
                let rot_elements = compute_rotational_elements(b1, self.epoch);
                let rot_matrix = rotation_matrix_inertial_to_fixed(rot_elements);
                let v_fixed = rot_matrix * self.vector;
                Ok(Coordinate::new(v_fixed, target_frame, self.epoch))
            },
            
            // ICRS -> BodyInertial
            (Frame::ICRS, Frame::BodyInertial(b)) => {
                let r_body = get_body_barycentric(b, self.epoch);
                Ok(Coordinate::new(self.vector - r_body, target_frame, self.epoch))
            },
            // BodyInertial -> ICRS
            (Frame::BodyInertial(b), Frame::ICRS) => {
                let r_body = get_body_barycentric(b, self.epoch);
                Ok(Coordinate::new(self.vector + r_body, target_frame, self.epoch))
            },

            // Recursive / Chained
            // BodyFixed -> ICRS (via BodyInertial)
            (Frame::BodyFixed(b), Frame::ICRS) => {
                let inertial = self.transform_to(Frame::BodyInertial(b))?;
                inertial.transform_to(Frame::ICRS)
            },
            // ICRS -> BodyFixed (via BodyInertial)
            (Frame::ICRS, Frame::BodyFixed(b)) => {
                let inertial = self.transform_to(Frame::BodyInertial(b))?;
                inertial.transform_to(Frame::BodyFixed(b))
            },

            // GCRS handling (alias for EarthInertial essentially)
            (Frame::GCRS, Frame::ICRS) => {
                 Coordinate::new(self.vector, Frame::BodyInertial(crate::bodies::EARTH), self.epoch)
                    .transform_to(Frame::ICRS)
                    .map(|mut c| { c.frame = Frame::ICRS; c })
            },
            (Frame::ICRS, Frame::GCRS) => {
                 Coordinate::new(self.vector, Frame::ICRS, self.epoch)
                    .transform_to(Frame::BodyInertial(crate::bodies::EARTH))
                    .map(|mut c| { c.frame = Frame::GCRS; c })
            },
            
            // GCRS <-> ITRS (BodyFixed(Earth))
            (Frame::GCRS, Frame::BodyFixed(b)) if b.name == "Earth" => {
                 Coordinate::new(self.vector, Frame::BodyInertial(b), self.epoch)
                    .transform_to(Frame::BodyFixed(b))
            },
             (Frame::BodyFixed(b), Frame::GCRS) if b.name == "Earth" => {
                 Coordinate::new(self.vector, Frame::BodyFixed(b), self.epoch)
                    .transform_to(Frame::BodyInertial(b))
                    .map(|mut c| { c.frame = Frame::GCRS; c })
            },

            _ => Err(format!("Transformation from {:?} to {:?} not implemented", self.frame, target_frame)),
        }
    }
}

fn rotation_matrix_inertial_to_fixed(elements: RotationalElements) -> Matrix3<f64> {
    let alpha_rad = elements.alpha_deg.to_radians();
    let delta_rad = elements.delta_deg.to_radians();
    let w_rad = elements.w_deg.to_radians();

    // W is rotation about Z (prime meridian)
    // alpha, delta define the pole.
    // Standard IAU:
    // M = R3(W) * R1(pi/2 - delta) * R3(pi/2 + alpha)
    
    let axis_z = Vector3::z_axis();
    let axis_x = Vector3::x_axis();

    let r3_w = Rotation3::from_axis_angle(&axis_z, w_rad);
    let r1_delta = Rotation3::from_axis_angle(&axis_x, FRAC_PI_2 - delta_rad);
    let r3_alpha = Rotation3::from_axis_angle(&axis_z, FRAC_PI_2 + alpha_rad);

    (r3_w * r1_delta * r3_alpha).into_inner()
}

// Placeholder for ephemeris
fn get_body_barycentric(_body: Body, _epoch: f64) -> Vector3<f64> {
    // TODO: Implement VSOP87 or similar. For now, zero vector.
    Vector3::zeros()
}

#[cfg(test)]
mod tests;
