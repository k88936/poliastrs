pub mod cr3bp;
pub mod flyby;
pub mod lagrange;
pub mod soi;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThreeBodyError {
    InvalidMassParameter,
    InvalidMassRatio,
    InvalidFlybyInput,
    NonConvergentSolver,
}
