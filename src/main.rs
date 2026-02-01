mod cli;
mod core;
mod report;
mod simd;

fn main() -> anyhow::Result<()> {
    cli::run::entry()
}
