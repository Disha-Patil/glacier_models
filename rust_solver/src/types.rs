//! Core data types and constants for the glacier solver

use ndarray::Array2;

/// Physical and numerical constants
pub struct Constants {
    /// Grid spacing in meters
    pub dx: f64,
    pub dy: f64,

    /// Time stepping
    pub dt: f64,
    pub tmax: f64,

    /// Physical parameters
    /// Glen's flow law parameter (Pa^-3 year^-1)
    pub c: f64,
    /// Mass balance gradient (year^-1 per meter elevation)
    pub mb_grad: f64,

    /// Derived constants
    pub inv_2_dx_dx: f64,
    pub c_inv_2_dx_sq: f64,
}

impl Constants {
    pub fn new(dx: f64, dy: f64, dt: f64, tmax: f64, c: f64, mb_grad: f64) -> Self {
        let inv_2_dx_dx = 1.0 / (2.0 * dx * dx);
        let inv_2_dx_sq = 1.0 / (2.0 * 2.0 * dx * dx);
        let c_inv_2_dx_sq = c * inv_2_dx_sq;

        Constants {
            dx,
            dy,
            dt,
            tmax,
            c,
            mb_grad,
            inv_2_dx_dx,
            c_inv_2_dx_sq,
        }
    }

    /// Default constants matching C++ implementation
    pub fn default_himalaya() -> Self {
        Constants::new(
            100.0,          // dx
            100.0,          // dy
            0.01,           // dt (years)
            10000.0,        // tmax (years)
            0.0000111,      // C (Glen's flow law)
            0.006,          // mb_grad
        )
    }
}

/// Solver state and statistics
#[derive(Debug, Clone)]
pub struct SolverStats {
    pub iteration: usize,
    pub time_years: f64,
    pub total_ice: f64,
    pub ice_gain: f64,
    pub ice_melt: f64,
    pub ice_b_melt: f64,
    pub area_fraction: f64,
    pub steady_slope: f64,
}

impl Default for SolverStats {
    fn default() -> Self {
        SolverStats {
            iteration: 0,
            time_years: 0.0,
            total_ice: 0.0,
            ice_gain: 0.0,
            ice_melt: 0.0,
            ice_b_melt: 0.0,
            area_fraction: 0.0,
            steady_slope: 0.0,
        }
    }
}

/// Neighbor information for finite difference stencil
#[derive(Debug, Clone, Copy)]
pub struct NeighborInfo {
    pub value: f64,  // surface elevation at neighbor
    pub mask: u8,    // 1 = central difference, 2 = forward/backward
}
