use nalgebra::Vector3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CartesianState {
    pub r_km: Vector3<f64>,
    pub v_km_s: Vector3<f64>,
}

impl CartesianState {
    pub fn new(r_km: Vector3<f64>, v_km_s: Vector3<f64>) -> Self {
        Self { r_km, v_km_s }
    }
}
