pub mod util;
pub mod coesa76;
// pub mod coesa62;
// pub mod jacchia77;

pub use coesa76::COESA76;

pub trait Atmosphere {
    fn density(&self, alt_km: f64) -> f64;
}
