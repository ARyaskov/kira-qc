use crate::core::engine::RunOutput;
use crate::core::model::Mode;
use crate::report::html;
use anyhow::{Context, Result};
use std::fs;
use std::io::Write;
use std::path::Path;
use svg2pdf::usvg;
use svg2pdf::{ConversionOptions, PageOptions};
use zip::ZipWriter;
use zip::write::FileOptions;

#[derive(Clone, Copy, Debug)]
pub enum LatexMode {
    Summary,
    Supplement,
}

pub fn write(out_dir: &Path, output: &RunOutput, mode: LatexMode) -> Result<()> {
    let metrics = output.agg.finalize(&output.ctx);
    let latex_dir = out_dir.join("latex");
    let figures_dir = latex_dir.join("figures");
    let tables_dir = latex_dir.join("tables");
    fs::create_dir_all(&figures_dir)?;
    fs::create_dir_all(&tables_dir)?;

    write_basic_stats_table(&tables_dir, &metrics, &output.ctx.file_name)?;

    let mut figures: Vec<Figure> = Vec::new();
    match mode {
        LatexMode::Summary => {
            if output.ctx.mode == Mode::Short {
                figures.push(fig(
                    "per_base_quality",
                    "Per base sequence quality",
                    html::latex_svg_per_base_quality(&metrics)?,
                ));
                figures.push(fig(
                    "per_sequence_quality",
                    "Per sequence quality scores",
                    html::latex_svg_per_seq_quality(&metrics)?,
                ));
                figures.push(fig(
                    "duplication_levels",
                    "Sequence duplication levels",
                    html::latex_svg_duplication(&metrics)?,
                ));
                if metrics.statuses.adapter_content != crate::core::model::Status::Pass {
                    figures.push(fig(
                        "adapter_content",
                        "Adapter content",
                        html::latex_svg_adapter_content(&metrics)?,
                    ));
                }
            } else {
                figures.push(fig(
                    "sequence_length_distribution",
                    "Sequence length distribution",
                    html::latex_svg_length_dist(&metrics)?,
                ));
                figures.push(fig(
                    "per_sequence_quality",
                    "Per sequence quality scores",
                    html::latex_svg_per_seq_quality(&metrics)?,
                ));
                if metrics.statuses.adapter_content != crate::core::model::Status::Pass {
                    figures.push(fig(
                        "adapter_content",
                        "Adapter content",
                        html::latex_svg_adapter_content(&metrics)?,
                    ));
                }
            }
        }
        LatexMode::Supplement => {
            if output.ctx.mode == Mode::Short {
                figures.extend([
                    fig(
                        "per_base_quality",
                        "Per base sequence quality",
                        html::latex_svg_per_base_quality(&metrics)?,
                    ),
                    fig(
                        "per_sequence_quality",
                        "Per sequence quality scores",
                        html::latex_svg_per_seq_quality(&metrics)?,
                    ),
                    fig(
                        "per_base_content",
                        "Per base sequence content",
                        html::latex_svg_per_base_content(&metrics)?,
                    ),
                    fig(
                        "per_sequence_gc",
                        "Per sequence GC content",
                        html::latex_svg_per_seq_gc(&metrics)?,
                    ),
                    fig(
                        "per_base_n",
                        "Per base N content",
                        html::latex_svg_per_base_n(&metrics)?,
                    ),
                    fig(
                        "sequence_length_distribution",
                        "Sequence length distribution",
                        html::latex_svg_length_dist(&metrics)?,
                    ),
                    fig(
                        "duplication_levels",
                        "Sequence duplication levels",
                        html::latex_svg_duplication(&metrics)?,
                    ),
                    fig(
                        "overrepresented_sequences",
                        "Overrepresented sequences",
                        html::latex_svg_overrep(&metrics)?,
                    ),
                    fig(
                        "adapter_content",
                        "Adapter content",
                        html::latex_svg_adapter_content(&metrics)?,
                    ),
                ]);
                #[cfg(not(feature = "no-kmer"))]
                figures.push(fig(
                    "kmer_content",
                    "Kmer content",
                    html::latex_svg_kmer_content(&metrics)?,
                ));
            } else {
                figures.extend([
                    fig(
                        "sequence_length_distribution",
                        "Sequence length distribution",
                        html::latex_svg_length_dist(&metrics)?,
                    ),
                    fig(
                        "per_sequence_quality",
                        "Per sequence quality scores",
                        html::latex_svg_per_seq_quality(&metrics)?,
                    ),
                    fig(
                        "per_sequence_gc",
                        "Per sequence GC content",
                        html::latex_svg_per_seq_gc(&metrics)?,
                    ),
                    fig(
                        "per_sequence_n",
                        "Per sequence N content",
                        html::latex_svg_per_seq_n(&metrics)?,
                    ),
                    fig(
                        "adapter_content",
                        "Adapter content",
                        html::latex_svg_adapter_content(&metrics)?,
                    ),
                ]);
            }
        }
    }

    write_figures(&figures_dir, &figures)?;
    write_readme(&latex_dir)?;
    write_tex(
        &latex_dir,
        &output.ctx.file_name,
        output.ctx.mode,
        mode,
        &figures,
    )?;
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

fn write_basic_stats_table(
    tables_dir: &Path,
    metrics: &crate::core::metrics::FinalMetrics,
    file: &str,
) -> Result<()> {
    let mut out = String::new();
    out.push_str("\\begin{tabular}{ll}\n");
    out.push_str("\\toprule\n");
    out.push_str("Measure & Value \\\\\n");
    out.push_str("\\midrule\n");
    out.push_str(&format!("Filename & {} \\\\\n", escape_tex(file)));
    out.push_str(&format!(
        "File type & {} \\\\\n",
        escape_tex(&metrics.basic.file_type)
    ));
    out.push_str(&format!(
        "Encoding & {} \\\\\n",
        escape_tex(&metrics.basic.encoding)
    ));
    out.push_str(&format!(
        "Total Sequences & {} \\\\\n",
        fmt_int(metrics.basic.total_sequences)
    ));
    out.push_str(&format!(
        "Filtered Sequences & {} \\\\\n",
        fmt_int(metrics.basic.filtered_sequences)
    ));
    if metrics.basic.min_len == metrics.basic.max_len {
        out.push_str(&format!(
            "Sequence length & {} \\\\\n",
            metrics.basic.min_len
        ));
    } else {
        out.push_str(&format!(
            "Sequence length & {}-{} \\\\\n",
            metrics.basic.min_len, metrics.basic.max_len
        ));
    }
    out.push_str(&format!("%GC & {} \\\\\n", metrics.basic.gc_percent));
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
    out.push_str("\\textbf{Tool:} kira-qc\n");
    out.push_str("\\section*{Basic statistics}\n");
    out.push_str("\\input{tables/basic_statistics.tex}\n");

    match export_mode {
        LatexMode::Summary => {
            out.push_str("\\section*{Quality metrics}\n");
            for f in figures {
                out.push_str("\\begin{figure}[H]\n");
                out.push_str("\\centering\n");
                out.push_str(&format!(
                    "\\includegraphics[width=\\linewidth]{{figures/{}.pdf}}\n",
                    f.name
                ));
                out.push_str(&format!("\\caption{{{}}}\n", f.caption));
                out.push_str("\\end{figure}\n");
            }
        }
        LatexMode::Supplement => {
            for f in figures {
                out.push_str(&format!("\\section*{{{}}}\n", f.caption));
                out.push_str("\\begin{figure}[H]\n");
                out.push_str("\\centering\n");
                out.push_str(&format!(
                    "\\includegraphics[width=\\linewidth]{{figures/{}.pdf}}\n",
                    f.name
                ));
                out.push_str(&format!("\\caption{{{}}}\n", f.caption));
                out.push_str("\\end{figure}\n");
            }
        }
    }
    out.push_str("\\end{document}\n");

    let path = latex_dir.join("kira_qc.tex");
    fs::write(&path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
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


or upload generated zip-file `kira_qc_latex.zip` to Overleaf.com
"#;
    let path = latex_dir.join("README.tex.md");
    fs::write(&path, content).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn fmt_int(v: u64) -> String {
    let s = v.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    let len = s.len();
    for (i, ch) in s.chars().enumerate() {
        if i != 0 && (len - i) % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out
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

fn svg_to_pdf(svg: &str) -> Result<Vec<u8>> {
    let mut opt = usvg::Options::default();
    opt.fontdb_mut().load_system_fonts();
    let tree =
        usvg::Tree::from_str(svg, &opt).map_err(|e| anyhow::anyhow!("usvg parse failed: {e}"))?;
    let pdf = svg2pdf::to_pdf(&tree, ConversionOptions::default(), PageOptions::default())
        .map_err(|e| anyhow::anyhow!("svg2pdf conversion failed: {e}"))?;
    Ok(pdf)
}

fn write_latex_zip(latex_dir: &Path) -> Result<()> {
    let zip_path = latex_dir.join("kira_qc_latex.zip");
    let file = fs::File::create(&zip_path)
        .with_context(|| format!("failed to create {}", zip_path.display()))?;
    let mut zip = ZipWriter::new(file);
    let opts: FileOptions<'static, ()> =
        FileOptions::default().compression_method(zip::CompressionMethod::Deflated);

    zip.add_directory("latex/", opts)?;
    zip.add_directory("latex/figures/", opts)?;
    zip.add_directory("latex/tables/", opts)?;

    add_file_to_zip(&mut zip, opts, latex_dir, "kira_qc.tex")?;
    add_file_to_zip(&mut zip, opts, latex_dir, "README.tex.md")?;
    add_file_to_zip(
        &mut zip,
        opts,
        &latex_dir.join("tables"),
        "basic_statistics.tex",
    )?;

    let figures_dir = latex_dir.join("figures");
    let mut entries = Vec::new();
    for entry in fs::read_dir(&figures_dir)? {
        let entry = entry?;
        if let Some(name) = entry.file_name().to_str().map(|s| s.to_string()) {
            entries.push(name);
        }
    }
    entries.sort();
    for name in entries {
        add_file_to_zip(&mut zip, opts, &figures_dir, &name)?;
    }

    zip.finish()?;
    Ok(())
}

fn add_file_to_zip(
    zip: &mut ZipWriter<fs::File>,
    opts: FileOptions<'static, ()>,
    dir: &Path,
    name: &str,
) -> Result<()> {
    let path = dir.join(name);
    let data = fs::read(&path).with_context(|| format!("failed to read {}", path.display()))?;
    let zip_name = if dir.ends_with("figures") {
        format!("latex/figures/{}", name)
    } else if dir.ends_with("tables") {
        format!("latex/tables/{}", name)
    } else {
        format!("latex/{}", name)
    };
    zip.start_file(zip_name, opts)?;
    zip.write_all(&data)?;
    Ok(())
}
