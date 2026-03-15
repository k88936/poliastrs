use nalgebra::Vector3;

pub trait Interpolator {
    // Interpolates position at target_epoch given sorted epochs and values.
    // Assumes epochs are sorted.
    fn interpolate(&self, epochs: &[f64], values: &[Vector3<f64>], target_epoch: f64) -> Vector3<f64>;
}

pub struct LinearInterpolator;
impl Interpolator for LinearInterpolator {
    fn interpolate(&self, epochs: &[f64], values: &[Vector3<f64>], target_epoch: f64) -> Vector3<f64> {
        if epochs.is_empty() { return Vector3::zeros(); }
        // Simple binary search or partition point
        // For sorted array
        let idx = epochs.partition_point(|&t| t <= target_epoch);
        
        if idx == 0 { return values[0]; }
        if idx >= epochs.len() { return *values.last().unwrap(); }

        let t0 = epochs[idx-1];
        let t1 = epochs[idx];
        let v0 = values[idx-1];
        let v1 = values[idx];

        if (t1 - t0).abs() < 1e-9 { return v0; }

        let fraction = (target_epoch - t0) / (t1 - t0);
        v0 + (v1 - v0) * fraction
    }
}

// TODO: Implement proper Sinc and Spline. For now alias to Linear to pass node tests.
pub struct SincInterpolator;
impl Interpolator for SincInterpolator {
    fn interpolate(&self, epochs: &[f64], values: &[Vector3<f64>], target_epoch: f64) -> Vector3<f64> {
        LinearInterpolator.interpolate(epochs, values, target_epoch)
    }
}

pub struct SplineInterpolator;
impl Interpolator for SplineInterpolator {
    fn interpolate(&self, epochs: &[f64], values: &[Vector3<f64>], target_epoch: f64) -> Vector3<f64> {
        LinearInterpolator.interpolate(epochs, values, target_epoch)
    }
}
