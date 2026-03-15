use nalgebra::Vector3;
use crate::bodies::Body;
use crate::frames::Plane;
use crate::twobody::orbit::Orbit;

pub mod interpolator;
use interpolator::{Interpolator, LinearInterpolator};

pub struct Ephem {
    pub epochs: Vec<f64>,
    pub coordinates: Vec<Vector3<f64>>,
    pub velocities: Option<Vec<Vector3<f64>>>,
    pub plane: Plane,
}

impl Ephem {
    pub fn new(epochs: Vec<f64>, coordinates: Vec<Vector3<f64>>, plane: Plane) -> Self {
        Ephem {
            epochs,
            coordinates,
            velocities: None,
            plane,
        }
    }

    pub fn from_orbit(orbit: Orbit, epochs: Vec<f64>, plane: Plane) -> Self {
        let mut coords = Vec::with_capacity(epochs.len());
        let mut vels = Vec::with_capacity(epochs.len());
        
        for &epoch in &epochs {
            let dt = epoch - orbit.epoch_tdb_seconds;
            let new_orbit = orbit.propagate_seconds(dt).unwrap();
            coords.push(new_orbit.state.r_km);
            vels.push(new_orbit.state.v_km_s);
        }
        
        Ephem {
            epochs,
            coordinates: coords,
            velocities: Some(vels),
            plane,
        }
    }

    pub fn from_body(_body: Body, epochs: Vec<f64>, plane: Plane) -> Self {
        // STUB: Return zeros for now as we lack VSOP87
        let coords = epochs.iter().map(|_| Vector3::zeros()).collect();
        Ephem {
            epochs,
            coordinates: coords,
            velocities: None,
            plane,
        }
    }

    pub fn sample<I: Interpolator>(&self, target_epochs: Option<Vec<f64>>, interpolator: I) -> Vec<Vector3<f64>> {
        let targets = target_epochs.unwrap_or_else(|| self.epochs.clone());
        targets.iter().map(|&t| interpolator.interpolate(&self.epochs, &self.coordinates, t)).collect()
    }

    pub fn rv(&self, epoch: Option<f64>) -> (Vec<Vector3<f64>>, Vec<Vector3<f64>>) {
        let targets = if let Some(t) = epoch { vec![t] } else { self.epochs.clone() };
        
        let interp = LinearInterpolator; 
        
        let r = targets.iter().map(|&t| interp.interpolate(&self.epochs, &self.coordinates, t)).collect();
        
        let v = if let Some(ref vels) = self.velocities {
            targets.iter().map(|&t| interp.interpolate(&self.epochs, vels, t)).collect()
        } else {
            targets.iter().map(|_| Vector3::zeros()).collect()
        };
        
        (r, v)
    }
}

#[cfg(test)]
mod tests;
