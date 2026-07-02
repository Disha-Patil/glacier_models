//! File I/O operations for the glacier solver

use ndarray::Array2;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// Read bedrock topography file
///
/// Expected format: space-separated values
/// i j elevation
pub fn read_bed_file(path: &Path) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut max_i = 0usize;
    let mut max_j = 0usize;
    let mut data = Vec::new();

    // First pass: find dimensions
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3 {
            let i: usize = parts[0].parse()?;
            let j: usize = parts[1].parse()?;
            let val: f64 = parts[2].parse()?;

            if i > max_i {
                max_i = i;
            }
            if j > max_j {
                max_j = j;
            }
            data.push((i, j, val));
        }
    }

    // Create matrix and fill values
    let mut matrix = Array2::zeros((max_i + 1, max_j + 1));
    for (i, j, val) in data {
        matrix[[i, j]] = val;
    }

    Ok(matrix)
}

/// Read a general field file (ELA, etc.)
pub fn read_field_file(path: &Path) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    read_bed_file(path)
}

/// Read ice mask file
pub fn read_ice_mask_file(path: &Path) -> Result<Array2<i32>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut max_i = 0usize;
    let mut max_j = 0usize;
    let mut data = Vec::new();

    // First pass: find dimensions
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3 {
            let i: usize = parts[0].parse()?;
            let j: usize = parts[1].parse()?;
            let val: i32 = parts[2].parse::<f64>()? as i32;

            if i > max_i {
                max_i = i;
            }
            if j > max_j {
                max_j = j;
            }
            data.push((i, j, val));
        }
    }

    // Create matrix and fill values
    let mut matrix = Array2::zeros((max_i + 1, max_j + 1));
    for (i, j, val) in data {
        matrix[[i, j]] = val;
    }

    Ok(matrix)
}

/// Read initial ice thickness file
pub fn read_h_init_file(path: &Path) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => {
            // If file doesn't exist, return zero matrix (will be sized later)
            return Ok(Array2::zeros((100, 100)));
        }
    };

    let reader = BufReader::new(file);
    let mut max_i = 0usize;
    let mut max_j = 0usize;
    let mut data = Vec::new();

    // First pass: find dimensions
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3 {
            let i: usize = parts[0].parse()?;
            let j: usize = parts[1].parse()?;
            let val: f64 = parts[2].parse()?;

            if i > max_i {
                max_i = i;
            }
            if j > max_j {
                max_j = j;
            }
            data.push((i, j, val));
        }
    }

    // Create matrix and fill values
    let mut matrix = Array2::zeros((max_i + 1, max_j + 1));
    for (i, j, val) in data {
        matrix[[i, j]] = val;
    }

    Ok(matrix)
}

/// Write final ice thickness and mass balance
pub fn write_h_steady(path: &Path, h: &Array2<f64>, mb: &Array2<f64>) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    let (xmax, ymax) = h.dim();

    for j in (0..ymax).rev() {
        for i in 0..xmax {
            writeln!(file, "{} {} {:.12e} {:.12e}", i, j, h[[i, j]], mb[[i, j]])?;
        }
    }

    Ok(())
}

/// Write steady state slope information
pub fn write_steady_slope(
    path: &Path,
    slope: f64,
    areaf: f64,
    tau: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(path)?;
    writeln!(file, "{:.12e} {:.6} {:.1}", slope, areaf, tau)?;
    Ok(())
}
