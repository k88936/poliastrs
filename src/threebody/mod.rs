pub mod cr3bp;
pub mod lagrange;
pub mod soi;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThreeBodyError {
    InvalidMassParameter,
    InvalidMassRatio,
    NonConvergentSolver,
}
