//! Matrix assembly functions for the glacier solver

use crate::types::{Constants, NeighborInfo};
use ndarray::{s, Array1, Array2};
use nalgebra_sparse::CsrMatrix;
use std::collections::BTreeMap;

/// Compute diffusion coefficient matrix
///
/// Based on the shallow ice approximation:
/// d = C * (1/2dx^2) * h^5 * (∇s)^2
///
pub fn fill_d(
    h: &Array2<f64>,
    bed: &Array2<f64>,
    m: &Array2<i32>,
    constants: &Constants,
) -> Array2<f64> {
    let (xmax, ymax) = h.dim();
    let mut d = Array2::zeros((xmax, ymax));
    let c_inv = constants.c_inv_2_dx_sq;

    // Extract neighbor information
    let mut neigh = vec![
        vec![NeighborInfo { value: 0.0, mask: 0 }; ymax]; // i+1
        4
    ];
    let mut ng = [0, 0, 0, 0];

    for i in 1..xmax - 1 {
        for j in 1..ymax - 1 {
            // i+1 neighbor
            if m[[i + 1, j]] != 0 {
                neigh[0].push(NeighborInfo {
                    value: bed[[i + 1, j]] + h[[i + 1, j]],
                    mask: 1,
                });
            } else {
                neigh[0].push(NeighborInfo {
                    value: bed[[i, j]] + h[[i, j]],
                    mask: 2,
                });
            }
            ng[0] += 1;

            // j+1 neighbor
            if m[[i, j + 1]] != 0 {
                neigh[1].push(NeighborInfo {
                    value: bed[[i, j + 1]] + h[[i, j + 1]],
                    mask: 1,
                });
            } else {
                neigh[1].push(NeighborInfo {
                    value: bed[[i, j]] + h[[i, j]],
                    mask: 2,
                });
            }
            ng[1] += 1;

            // i-1 neighbor
            if m[[i - 1, j]] != 0 {
                neigh[2].push(NeighborInfo {
                    value: bed[[i - 1, j]] + h[[i - 1, j]],
                    mask: 1,
                });
            } else {
                neigh[2].push(NeighborInfo {
                    value: bed[[i, j]] + h[[i, j]],
                    mask: 2,
                });
            }
            ng[2] += 1;

            // j-1 neighbor
            if m[[i, j - 1]] != 0 {
                neigh[3].push(NeighborInfo {
                    value: bed[[i, j - 1]] + h[[i, j - 1]],
                    mask: 1,
                });
            } else {
                neigh[3].push(NeighborInfo {
                    value: bed[[i, j]] + h[[i, j]],
                    mask: 2,
                });
            }
            ng[3] += 1;
        }
    }

    // Compute diffusion coefficient
    let mut ngg = 0;
    for i in 1..xmax - 1 {
        for j in 1..ymax - 1 {
            let h_power = h[[i, j]].powi(5);
            let grad_x_sq =
                (neigh[0][ngg].mask as f64 * neigh[2][ngg].mask as f64 * (neigh[0][ngg].value - neigh[2][ngg].value)).powi(2);
            let grad_y_sq =
                (neigh[1][ngg].mask as f64 * neigh[3][ngg].mask as f64 * (neigh[1][ngg].value - neigh[3][ngg].value)).powi(2);

            d[[i, j]] = c_inv * h_power * (grad_x_sq + grad_y_sq);
            ngg += 1;
        }
    }

    d
}

/// Compute mass balance
///
/// Linear mass balance model:
/// mb = mb_grad * (surface_elevation - ELA)
///
pub fn fill_mb(
    bed: &Array2<f64>,
    h: &Array2<f64>,
    ela: &Array2<f64>,
    domain: &[(usize, usize)],
    mb_grad: f64,
) -> Array2<f64> {
    let (xmax, ymax) = bed.dim();
    let mut mb = Array2::zeros((xmax, ymax));

    for &(i, j) in domain {
        let surface = bed[[i, j]] + h[[i, j]];
        let q = mb_grad * (surface - ela[[i, j]]);
        // Limit accumulation to 1.0
        mb[[i, j]] = if q >= 1.0 { 1.0 } else { q };
    }

    mb
}

/// Calculate area fraction and total area
pub fn calc_areaf(
    h: &Array2<f64>,
    ice_mask: &Array2<i32>,
) -> (f64, f64) {
    let mut nu = 0usize;
    let mut den = 0usize;

    for ((i, j), &h_val) in h.indexed_iter() {
        if h_val > 0.0 && ice_mask[[i, j]] > 0 {
            nu += 1;
        }
        if ice_mask[[i, j]] > 0 {
            den += 1;
        }
    }

    let areaf = nu as f64 / den as f64;
    (areaf, nu as f64)
}

/// Fill RHS vector b for linear system Ax = b
pub fn fill_b(
    h: &Array2<f64>,
    bed: &Array2<f64>,
    d: &Array2<f64>,
    mb: &Array2<f64>,
    m: &Array2<i32>,
    domain: &[(usize, usize)],
    dt: f64,
    dt_inv_2_dx_dx: f64,
) -> Array1<f64> {
    let mut b = Array1::zeros(domain.len());

    for (k, &(i, j)) in domain.iter().enumerate() {
        let alpha = d[[i + 1, j]] + d[[i, j]];
        let beta = d[[i, j]] + d[[i - 1, j]];
        let gamma = d[[i, j + 1]] + d[[i, j]];
        let eps = d[[i, j]] + d[[i, j - 1]];

        let flux = dt_inv_2_dx_dx
            * (alpha * (bed[[i + 1, j]] - bed[[i, j]]) * m[[i + 1, j]] as f64
                - beta * (bed[[i, j]] - bed[[i - 1, j]]) * m[[i - 1, j]] as f64
                + gamma * (bed[[i, j + 1]] - bed[[i, j]]) * m[[i, j + 1]] as f64
                - eps * (bed[[i, j]] - bed[[i, j - 1]]) * m[[i, j - 1]] as f64
                + (alpha * (h[[i + 1, j]] - h[[i, j]]) * m[[i + 1, j]] as f64
                    - beta * (h[[i, j]] - h[[i - 1, j]]) * m[[i - 1, j]] as f64
                    + gamma * (h[[i, j + 1]] - h[[i, j]]) * m[[i, j + 1]] as f64
                    - eps * (h[[i, j]] - h[[i, j - 1]]) * m[[i, j - 1]] as f64)
                    * 0.5);

        b[k] = h[[i, j]]
            + dt * mb[[i, j]]
            + dt_inv_2_dx_dx
                * (alpha * (bed[[i + 1, j]] - bed[[i, j]]) * m[[i + 1, j]] as f64
                    - beta * (bed[[i, j]] - bed[[i - 1, j]]) * m[[i - 1, j]] as f64
                    + gamma * (bed[[i, j + 1]] - bed[[i, j]]) * m[[i, j + 1]] as f64
                    - eps * (bed[[i, j]] - bed[[i, j - 1]]) * m[[i, j - 1]] as f64
                    + (alpha * (h[[i + 1, j]] - h[[i, j]]) * m[[i + 1, j]] as f64
                        - beta * (h[[i, j]] - h[[i - 1, j]]) * m[[i - 1, j]] as f64
                        + gamma * (h[[i, j + 1]] - h[[i, j]]) * m[[i, j + 1]] as f64
                        - eps * (h[[i, j]] - h[[i, j - 1]]) * m[[i, j - 1]] as f64)
                        * 0.5);
    }

    b
}

/// Build sparse system matrix A for implicit time stepping
pub fn fill_A(
    d: &Array2<f64>,
    m: &Array2<i32>,
    m_idx: &Array2<i32>,
    domain: &[(usize, usize)],
    dt_inv_2_dx_dx: f64,
) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
    // Build COO format (row, col, value) for sparse matrix
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut values = Vec::new();

    for (l, &(i, j)) in domain.iter().enumerate() {
        let alpha = d[[i + 1, j]] + d[[i, j]];
        let beta = d[[i, j]] + d[[i - 1, j]];
        let gamma = d[[i, j + 1]] + d[[i, j]];
        let eps = d[[i, j]] + d[[i, j - 1]];

        let c1 = 1.0 + (dt_inv_2_dx_dx * (alpha * m[[i + 1, j]] as f64
            + beta * m[[i - 1, j]] as f64
            + gamma * m[[i, j + 1]] as f64
            + eps * m[[i, j - 1]] as f64)
            * 0.5);

        // Diagonal element
        rows.push(l);
        cols.push(l);
        values.push(c1);

        // Off-diagonal elements
        if m[[i + 1, j]] != 0 {
            rows.push(l);
            cols.push(m_idx[[i + 1, j]] as usize);
            values.push(-dt_inv_2_dx_dx * alpha * m[[i + 1, j]] as f64 * 0.5);
        }
        if m[[i - 1, j]] != 0 {
            rows.push(l);
            cols.push(m_idx[[i - 1, j]] as usize);
            values.push(-dt_inv_2_dx_dx * beta * m[[i - 1, j]] as f64 * 0.5);
        }
        if m[[i, j + 1]] != 0 {
            rows.push(l);
            cols.push(m_idx[[i, j + 1]] as usize);
            values.push(-dt_inv_2_dx_dx * gamma * m[[i, j + 1]] as f64 * 0.5);
        }
        if m[[i, j - 1]] != 0 {
            rows.push(l);
            cols.push(m_idx[[i, j - 1]] as usize);
            values.push(-dt_inv_2_dx_dx * eps * m[[i, j - 1]] as f64 * 0.5);
        }
    }

    (rows, cols, values)
}
