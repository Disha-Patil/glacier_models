//! Main glacier solver implementation

use crate::io;
use crate::matrices::{calc_areaf, fill_A, fill_b, fill_d, fill_mb};
use crate::types::{Constants, SolverStats};
use ndarray::{s, Array1, Array2, Array2 as Matrix};
use std::collections::VecDeque;

/// Main glacier solver struct
pub struct GlacierSolver {
    // Field matrices
    h: Array2<f64>,        // ice thickness
    bed: Array2<f64>,      // bedrock topography
    bed0: Array2<f64>,     // bedrock without smoothing
    ela: Array2<f64>,      // equilibrium line altitude
    mb: Array2<f64>,       // mass balance
    d: Array2<f64>,        // diffusion coefficient
    ice_mask: Array2<i32>, // glacier domain mask
    m: Array2<i32>,        // domain mask (same as ice_mask for now)
    m_idx: Array2<i32>,    // sequential index for domain points

    // Domain information
    domain: Vec<(usize, usize)>,
    xmax: usize,
    ymax: usize,
    n0: usize, // number of domain points

    // Constants
    constants: Constants,

    // Statistics
    stats: SolverStats,
}

impl GlacierSolver {
    /// Create a new solver instance
    pub fn new(
        bed: Array2<f64>,
        bed0: Array2<f64>,
        ela: Array2<f64>,
        ice_mask: Array2<i32>,
        h_init: Array2<f64>,
    ) -> Self {
        let (xmax, ymax) = bed.dim();
        let constants = Constants::default_himalaya();

        // Create mask from bedrock (non-zero = glacier domain)
        let mut m = Array2::zeros((xmax, ymax));
        for ((i, j), &bed_val) in bed.indexed_iter() {
            if bed_val > 0.0 {
                m[[i, j]] = 1;
            }
        }

        // Extract domain points
        let mut domain = Vec::new();
        for i in 0..xmax {
            for j in 0..ymax {
                if m[[i, j]] != 0 {
                    domain.push((i, j));
                }
            }
        }
        let n0 = domain.len();

        // Create sequential index map
        let mut m_idx = Array2::zeros((xmax, ymax));
        for (idx, &(i, j)) in domain.iter().enumerate() {
            m_idx[[i, j]] = idx as i32;
        }

        println!("  Domain: {} x {} grid, {} active points", xmax, ymax, n0);

        GlacierSolver {
            h: h_init,
            bed,
            bed0,
            ela,
            mb: Array2::zeros((xmax, ymax)),
            d: Array2::zeros((xmax, ymax)),
            ice_mask,
            m,
            m_idx,
            domain,
            xmax,
            ymax,
            n0,
            constants,
            stats: SolverStats::default(),
        }
    }

    /// Return domain points
    pub fn domain_points(&self) -> &[(usize, usize)] {
        &self.domain
    }

    /// Main solver loop
    pub fn solve(&mut self) -> Result<(Array2<f64>, Array2<f64>, f64, f64, f64), Box<dyn std::error::Error>> {
        let mut time = 0.0;
        let mut count = 0usize;
        let mut ice_count = 0usize;
        let mut ice_arr = [0.0; 2];
        let mut areaf_arr = [0.0; 2];

        let mut prev_totalice = 0.0;
        let mut prev_gain = 0.0;
        let mut prev_melt = 0.0;
        let mut prev_b_melt = 0.0;

        let mut tau = 0.0;
        let mut steady_slope = 0.0;
        let mut areaf = 0.0;
        let mut area = 0.0;

        let dt = self.constants.dt;
        let tmax = self.constants.tmax;
        let dt_inv_2_dx_dx = dt * self.constants.inv_2_dx_dx;

        let mut ice_gain = 0.0;
        let mut ice_melt = 0.0;
        let mut ice_b_melt = 0.0;

        println!("\n⏳ Starting time loop...");
        println!("   dt = {} years, tmax = {} years", dt, tmax);

        while time <= tmax {
            count += 1;
            ice_count += 1;

            // Calculate area fraction
            (areaf, area) = calc_areaf(&self.h, &self.ice_mask);

            if areaf >= 0.8 && tau == 0.0 {
                tau = ((time / 25.0).floor() + 1.0) * 25.0;
            }

            // Break conditions
            if time == 5000.0 {
                tau = 1500.0;
                println!("   ⏸️  Break: reached tmax intermediate check (5000 years)");
                break;
            }

            if time > 200.0
                && areaf_arr[0] == areaf_arr[1]
                && areaf < 0.8
                && steady_slope < 0.000001
            {
                tau = time;
                println!("   ⏸️  Break: never reached 80% coverage");
                break;
            }

            // Update mass balance
            self.mb = fill_mb(&self.bed, &self.h, &self.ela, &self.domain, self.constants.mb_grad);

            // Update diffusion
            self.d = fill_d(&self.h, &self.bed, &self.m, &self.constants);

            // Fill RHS vector and system matrix
            let b = fill_b(
                &self.h,
                &self.bed,
                &self.d,
                &self.mb,
                &self.m,
                &self.domain,
                dt,
                dt_inv_2_dx_dx,
            );

            let (rows, cols, values) =
                fill_A(&self.d, &self.m, &self.m_idx, &self.domain, dt_inv_2_dx_dx);

            // Solve linear system using conjugate gradient
            let x = self.solve_cg(&b, &rows, &cols, &values)?;

            // Update h matrix and compute statistics
            let mut totalice = 0.0;
            let mut net_mb = 0.0;
            let mut inst_area = 0usize;

            ice_gain = 0.0;
            ice_melt = 0.0;
            ice_b_melt = 0.0;

            for (p, &(i, j)) in self.domain.iter().enumerate() {
                if x[p] > 0.0 {
                    if self.mb[[i, j]] > 0.0 {
                        ice_gain += dt * self.mb[[i, j]];
                    } else {
                        ice_melt += dt * self.mb[[i, j]];
                    }
                } else {
                    ice_b_melt += -dt * self.mb[[i, j]] + x[p];
                }

                self.h[[i, j]] = x[p].max(0.0); // negative ice → 0

                totalice += self.h[[i, j]];
                if self.h[[i, j]] > 0.0 {
                    net_mb += self.mb[[i, j]];
                    inst_area += 1;
                }
            }

            // Check for steady state every 20000 iterations
            if count == 20000 {
                count = 0;
                areaf_arr[0] = areaf_arr[1];
                areaf_arr[1] = areaf;

                steady_slope = if inst_area > 0 {
                    (totalice - prev_totalice) / (inst_area as f64 * 200.0)
                } else {
                    0.0
                };

                if areaf_arr[0] == areaf_arr[1] && steady_slope > 0.8 {
                    tau = time;
                    println!("   ⏸️  Break: ice accumulation too fast");
                    break;
                }

                if areaf_arr[0] == areaf_arr[1] && steady_slope < 0.0001 {
                    tau = time;
                    println!("   ⏸️  Break: reached steady state");
                    break;
                }
            }

            // Log progress every 20000 ice count iterations
            if ice_count == 20000 {
                ice_count = 0;
                println!(
                    "   t={:.0} yrs | slope={:.2e} | areaf={:.4} | totalice={:.2e}",
                    time, steady_slope, areaf, totalice
                );
            }

            prev_gain = ice_gain;
            prev_melt = ice_melt;
            prev_b_melt = ice_b_melt;
            prev_totalice = totalice;

            time += dt;
        }

        println!(
            "\n✓ Steady state reached: slope={:.2e}, areaf={:.4}, tau={:.1} years",
            steady_slope, areaf, tau
        );

        Ok((self.h.clone(), self.mb.clone(), steady_slope, areaf, tau))
    }

    /// Solve linear system using Conjugate Gradient method
    fn solve_cg(
        &self,
        b: &Array1<f64>,
        rows: &[usize],
        cols: &[usize],
        values: &[f64],
    ) -> Result<Array1<f64>, Box<dyn std::error::Error>> {
        // Simple CG implementation (naive version)
        // For production, consider using a dedicated sparse solver library
        let n = b.len();
        let mut x = Array1::zeros(n);

        // Build sparse matrix function
        let apply_A = |v: &Array1<f64>| -> Array1<f64> {
            let mut result = Array1::zeros(n);
            for (idx, (r, c)) in rows.iter().zip(cols.iter()).enumerate() {
                result[*r] += values[idx] * v[*c];
            }
            result
        };

        // CG iterations
        let mut r = b - &apply_A(&x);
        let mut p = r.clone();
        let mut rsold = (&r * &r).sum();
        let tolerance = 1e-9;
        let max_iters = 1000;

        for iter in 0..max_iters {
            let ap = apply_A(&p);
            let pap = (&p * &ap).sum();

            if pap.abs() < 1e-14 {
                break;
            }

            let alpha = rsold / pap;
            x = &x + alpha * &p;
            r = &r - alpha * &ap;
            let rsnew = (&r * &r).sum();

            if rsnew.sqrt() < tolerance {
                break;
            }

            let beta = rsnew / rsold;
            p = &r + beta * &p;
            rsold = rsnew;
        }

        Ok(x)
    }
}
