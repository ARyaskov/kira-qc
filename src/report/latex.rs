use crate::core::metrics::{BasicStats, FinalMetrics};
use crate::core::model::{FinalizeContext, Mode, Status};
use crate::report::zip::{deterministic_time, list_files_recursive};
use crate::report::{ModuleId, active_modules, html};
use anyhow::{Context, Result};
use std::fs;
use std::io::Write;
use std::path::Path;
use std::sync::OnceLock;
use svg2pdf::usvg;
use svg2pdf::{ConversionOptions, PageOptions};
use zip::ZipWriter;
use zip::write::FileOptions;

#[derive(Clone, Copy, Debug)]
pub enum LatexMode {
    Summary,
    Supplement,
}

pub fn write(
    out_dir: &Path,
    metrics: &FinalMetrics,
    ctx: &FinalizeContext,
    mode: LatexMode,
) -> Result<()> {
    let latex_dir = out_dir.join("latex");
    let figures_dir = latex_dir.join("figures");
    let tables_dir = latex_dir.join("tables");
    fs::create_dir_all(&figures_dir)?;
    fs::create_dir_all(&tables_dir)?;

    let active = active_modules(metrics, ctx);
    let basic = active.contains(&ModuleId::Basic);
    if basic {
        write_basic_stats_table(&tables_dir, metrics, &ctx.file_name)?;
    }

    let figures = collect_figures(metrics, ctx, mode, &active)?;
    write_figures(&figures_dir, &figures)?;
    write_readme(&latex_dir)?;
    write_tex(&latex_dir, &ctx.file_name, ctx.mode, mode, &figures, basic)?;
    write_latex_zip(&latex_dir)?;
    Ok(())
}

struct Figure {
    name: &'static str,
    caption: &'static str,
    svg: String,
}

fn fig(name: &'static str, caption: &'static str, svg: String) -> Figure {
    Figure { name, caption, svg }
}

fn collect_figures(
    metrics: &FinalMetrics,
    ctx: &FinalizeContext,
    mode: LatexMode,
    active: &[ModuleId],
) -> Result<Vec<Figure>> {
    use ModuleId::*;
    let selected = match (ctx.mode, mode) {
        (Mode::Short, LatexMode::Summary) => {
            &[PerBaseQual, PerSeqQual, Duplication, AdapterContent][..]
        }
        (Mode::Long, LatexMode::Summary) => &[LengthDist, PerSeqQual, AdapterContent][..],
        (_, LatexMode::Supplement) => active,
    };
    selected
        .iter()
        .copied()
        .filter(|module| *module != Basic && active.contains(module))
        .filter(|module| {
            !matches!(mode, LatexMode::Summary)
                || *module != AdapterContent
                || metrics.statuses.adapter_content != Status::Pass
        })
        .map(|module| {
            let (name, caption, render): (_, _, fn(&FinalMetrics) -> Result<String>) = match module
            {
                Basic => unreachable!("basic statistics are rendered as a table"),
                PerBaseQual => (
                    "per_base_quality",
                    "Per base sequence quality",
                    html::latex_svg_per_base_quality,
                ),
                PerSeqQual => (
                    "per_sequence_quality",
                    "Per sequence quality scores",
                    html::latex_svg_per_seq_quality,
                ),
                Duplication => (
                    "duplication_levels",
                    "Sequence duplication levels",
                    html::latex_svg_duplication,
                ),
                AdapterContent => (
                    "adapter_content",
                    "Adapter content",
                    html::latex_svg_adapter_content,
                ),
                LengthDist => (
                    "sequence_length_distribution",
                    "Sequence length distribution",
                    html::latex_svg_length_dist,
                ),
                PerTile => (
                    "per_tile_quality",
                    "Per tile sequence quality",
                    html::latex_svg_per_tile,
                ),
                PerBaseContent => (
                    "per_base_content",
                    "Per base sequence content",
                    html::latex_svg_per_base_content,
                ),
                PerSeqGc => (
                    "per_sequence_gc",
                    "Per sequence GC content",
                    html::latex_svg_per_seq_gc,
                ),
                PerBaseN => (
                    "per_base_n",
                    "Per base N content",
                    html::latex_svg_per_base_n,
                ),
                Overrepresented => (
                    "overrepresented_sequences",
                    "Overrepresented sequences",
                    html::latex_svg_overrep,
                ),
                KmerContent => ("kmer_content", "Kmer content", html::latex_svg_kmer_content),
                PerSeqN => (
                    "per_sequence_n",
                    "Per sequence N content",
                    html::latex_svg_per_seq_n,
                ),
            };
            Ok(fig(name, caption, render(metrics)?))
        })
        .collect()
}

fn write_figures(dir: &Path, figures: &[Figure]) -> Result<()> {
    for f in figures {
        let svg_path = dir.join(format!("{}.svg", f.name));
        fs::write(&svg_path, &f.svg)
            .with_context(|| format!("failed to write {}", svg_path.display()))?;
        let pdf =
            svg_to_pdf(&f.svg).with_context(|| format!("failed to convert {} to PDF", f.name))?;
        let pdf_path = dir.join(format!("{}.pdf", f.name));
        fs::write(&pdf_path, pdf)
            .with_context(|| format!("failed to write {}", pdf_path.display()))?;
    }
    Ok(())
}

fn write_basic_stats_table(tables_dir: &Path, metrics: &FinalMetrics, file: &str) -> Result<()> {
    let b = &metrics.basic;
    let mut out = String::new();
    out.push_str("\\begin{tabular}{ll}\n");
    out.push_str("\\toprule\n");
    out.push_str("Measure & Value \\\\\n");
    out.push_str("\\midrule\n");
    out.push_str(&format!("Filename & {} \\\\\n", escape_tex(file)));
    out.push_str(&format!("File type & {} \\\\\n", escape_tex(b.file_type)));
    out.push_str(&format!("Encoding & {} \\\\\n", escape_tex(b.encoding)));
    out.push_str(&format!(
        "Total Sequences & {} \\\\\n",
        escape_tex(&crate::report::fmt_int(b.total_sequences))
    ));
    out.push_str(&format!(
        "Total Bases & {} \\\\\n",
        escape_tex(&BasicStats::format_bases(b.total_bases))
    ));
    out.push_str(&format!(
        "Sequences flagged as poor quality & {} \\\\\n",
        escape_tex(&crate::report::fmt_int(b.poor_quality_sequences))
    ));
    out.push_str(&format!(
        "Sequence length & {} \\\\\n",
        escape_tex(&b.length_range())
    ));
    out.push_str(&format!("\\%GC & {} \\\\\n", b.gc_percent));
    out.push_str("\\bottomrule\n");
    out.push_str("\\end{tabular}\n");
    let path = tables_dir.join("basic_statistics.tex");
    fs::write(&path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn write_tex(
    latex_dir: &Path,
    file_name: &str,
    mode: Mode,
    export_mode: LatexMode,
    figures: &[Figure],
    basic: bool,
) -> Result<()> {
    let mut out = String::new();
    out.push_str("\\documentclass{article}\n");
    out.push_str("\\usepackage{graphicx}\n");
    out.push_str("\\usepackage{booktabs}\n");
    out.push_str("\\usepackage{caption}\n");
    out.push_str("\\usepackage{float}\n");
    out.push_str("\\usepackage{geometry}\n");
    out.push_str("\\geometry{margin=1in}\n");
    out.push_str("\\title{Quality Control Report}\n");
    out.push_str("\\author{kira-qc}\n");
    out.push_str("\\date{\\today}\n");
    out.push_str("\\begin{document}\n");
    out.push_str("\\maketitle\n");
    out.push_str("\\section*{Sample information}\n");
    out.push_str(&format!(
        "\\textbf{{Input:}} {}\\\\\n",
        escape_tex(file_name)
    ));
    out.push_str(&format!(
        "\\textbf{{Mode:}} {}\\\\\n",
        match mode {
            Mode::Short => "Short-read (Illumina)",
            Mode::Long => "Long-read (ONT / PacBio)",
        }
    ));
    out.push_str(&format!(
        "\\textbf{{Tool:}} kira-qc {}\n",
        escape_tex(crate::KIRA_QC_VERSION)
    ));
    if basic {
        out.push_str("\\section*{Basic statistics}\n");
        out.push_str("\\input{tables/basic_statistics.tex}\n");
    }

    match export_mode {
        LatexMode::Summary => {
            out.push_str("\\section*{Quality metrics}\n");
            for f in figures {
                push_figure(&mut out, f);
            }
        }
        LatexMode::Supplement => {
            for f in figures {
                out.push_str(&format!("\\section*{{{}}}\n", escape_tex(f.caption)));
                push_figure(&mut out, f);
            }
        }
    }
    out.push_str("\\end{document}\n");

    let path = latex_dir.join("kira_qc.tex");
    fs::write(&path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn push_figure(out: &mut String, f: &Figure) {
    out.push_str("\\begin{figure}[H]\n");
    out.push_str("\\centering\n");
    out.push_str(&format!(
        "\\includegraphics[width=\\linewidth]{{figures/{}.pdf}}\n",
        f.name
    ));
    out.push_str(&format!("\\caption{{{}}}\n", escape_tex(f.caption)));
    out.push_str("\\end{figure}\n");
}

fn write_readme(latex_dir: &Path) -> Result<()> {
    let content = r#"# LaTeX export (kira-qc)

This folder contains a self-contained LaTeX report (kira_qc.tex), SVG/PDF figures, and LaTeX tables.

PDF figures are generated directly by kira-qc. If you need to regenerate them from SVG, you can use ImageMagick:

```
magick figures/per_base_quality.svg figures/per_base_quality.pdf
```

## Compile the report

```
pdflatex kira_qc.tex
```

or upload the generated `kira_qc_latex.zip` to Overleaf.
"#;
    let path = latex_dir.join("README.md");
    fs::write(&path, content).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn escape_tex(s: &str) -> String {
    s.replace('\\', "\\textbackslash{}")
        .replace('&', "\\&")
        .replace('%', "\\%")
        .replace('$', "\\$")
        .replace('#', "\\#")
        .replace('_', "\\_")
        .replace('{', "\\{")
        .replace('}', "\\}")
        .replace('~', "\\textasciitilde{}")
        .replace('^', "\\textasciicircum{}")
}

/// System font loading takes hundreds of milliseconds; doing it per figure made
/// `--export-latex supplement` spend most of its time rescanning font
/// directories.
fn usvg_options() -> &'static usvg::Options<'static> {
    static OPTS: OnceLock<usvg::Options<'static>> = OnceLock::new();
    OPTS.get_or_init(|| {
        let mut opt = usvg::Options::default();
        opt.fontdb_mut().load_system_fonts();
        opt
    })
}

fn svg_to_pdf(svg: &str) -> Result<Vec<u8>> {
    let tree = usvg::Tree::from_str(svg, usvg_options())
        .map_err(|e| anyhow::anyhow!("usvg parse failed: {e}"))?;
    let pdf = svg2pdf::to_pdf(&tree, ConversionOptions::default(), PageOptions::default())
        .map_err(|e| anyhow::anyhow!("svg2pdf conversion failed: {e}"))?;
    Ok(pdf)
}

fn write_latex_zip(latex_dir: &Path) -> Result<()> {
    let zip_name = "kira_qc_latex.zip";
    let zip_path = latex_dir.join(zip_name);
    let file = fs::File::create(&zip_path)
        .with_context(|| format!("failed to create {}", zip_path.display()))?;
    let mut zip = ZipWriter::new(file);
    // Fixed timestamps: the main bundle is already reproducible and this one
    // used to embed the current clock.
    let opts: FileOptions<'static, ()> = FileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        .last_modified_time(deterministic_time());

    zip.add_directory("latex/", opts)?;
    zip.add_directory("latex/figures/", opts)?;
    zip.add_directory("latex/tables/", opts)?;

    let mut names = list_files_recursive(latex_dir, "")?;
    names.retain(|n| n != zip_name);
    names.sort();

    for name in names {
        let path = latex_dir.join(name.replace('/', std::path::MAIN_SEPARATOR_STR));
        let data = fs::read(&path).with_context(|| format!("failed to read {}", path.display()))?;
        zip.start_file(format!("latex/{}", name), opts)?;
        zip.write_all(&data)?;
    }

    zip.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tex_special_characters_are_escaped() {
        assert_eq!(escape_tex("a_b"), "a\\_b");
        assert_eq!(escape_tex("100%"), "100\\%");
        assert_eq!(escape_tex("a&b"), "a\\&b");
    }
}
