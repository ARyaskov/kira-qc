use crate::core::metrics::FinalMetrics;
use crate::core::model::FinalizeContext;
use crate::report::{active_modules, write_buffered};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub fn write(path: &Path, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    write_buffered(file, |writer| write_to(writer, metrics, ctx))
        .with_context(|| format!("failed to write {}", path.display()))
}

pub fn write_to(w: &mut dyn Write, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    for module in active_modules(metrics, ctx) {
        writeln!(
            w,
            "{}\t{}\t{}",
            module.status(metrics).as_str_upper(),
            module.title(),
            ctx.file_name
        )?;
    }
    Ok(())
}
