pub mod lambert;
pub mod vallado;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LambertError {
    // Izzo errors
    CollinearVectors,
    NoFeasibleSolution,
    ConvergenceFailed,
    DerivativeZero,
    TimeOfFlightNegative,
    
    // Vallado errors
    MultiRevolutionUnsupported,
    PhaseAngle180Deg,
    MaximumIterationsReached,
}

pub use lambert::izzo;

