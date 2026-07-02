mod solver;
mod matrices;
mod io;
mod types;

use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "Glacier Solver")]
#[command(about = "Shallow Ice Approximation (SIA) glacier model solver")]
struct Args {
    /// Input bed file (smoothbed.txt)
    #[arg(short, long, default_value = "smoothbed.txt")]
    bed_file: PathBuf,

    /// Input bed file without smoothing (bed.txt)
    #[arg(short, long, default_value = "bed.txt")]
    bed0_file: PathBuf,

    /// ELA file (ela.txt)
    #[arg(short, long, default_value = "ela.txt")]
    ela_file: PathBuf,

    /// Ice mask file (ice_mask)
    #[arg(short, long, default_value = "ice_mask")]
    ice_mask_file: PathBuf,

    /// Initial ice thickness file (h_init_file.txt)
    #[arg(short, long, default_value = "h_init_file.txt")]
    init_file: PathBuf,

    /// Output file for results (h_steady.txt)
    #[arg(short, long, default_value = "h_steady.txt")]
    output_file: PathBuf,

    /// Steady slope output file
    #[arg(short, long, default_value = "steady_slope")]
    slope_file: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    println!("🏔️  Glacier Solver (Rust Implementation)");
    println!("================================================");

    // Read input files
    println!("\n📂 Reading input files...");
    let bed = io::read_bed_file(&args.bed_file)?;
    let bed0 = io::read_bed_file(&args.bed0_file)?;
    let ela = io::read_field_file(&args.ela_file)?;
    let ice_mask = io::read_ice_mask_file(&args.ice_mask_file)?;
    let h_init = io::read_h_init_file(&args.init_file)?;

    println!("✓ Bed shape: {} x {}", bed.nrows(), bed.ncols());
    println!("✓ ELA shape: {} x {}", ela.nrows(), ela.ncols());

    // Initialize solver
    println!("\n🔧 Initializing solver...");
    let mut solver = solver::GlacierSolver::new(
        bed.clone(),
        bed0.clone(),
        ela.clone(),
        ice_mask.clone(),
        h_init,
    );

    println!("✓ Domain points: {}", solver.domain_points().len());

    // Run solver
    println!("\n⏱️  Running transient simulation...");
    let start = std::time::Instant::now();
    let (final_h, final_mb, slope, areaf, tau) = solver.solve()?;
    let elapsed = start.elapsed();

    println!("✓ Simulation complete in {:.2}s", elapsed.as_secs_f64());
    println!("  - Final steady slope: {:.12e}", slope);
    println!("  - Area fraction: {:.6}", areaf);
    println!("  - Time to steady: {:.1} years", tau);

    // Write output files
    println!("\n💾 Writing output files...");
    io::write_h_steady(&args.output_file, &final_h, &final_mb)?;
    io::write_steady_slope(&args.slope_file, slope, areaf, tau)?;

    println!("✓ Output written to {}", args.output_file.display());
    println!("✓ Slope written to {}", args.slope_file.display());

    println!("\n✅ Solver finished successfully!");
    Ok(())
}
