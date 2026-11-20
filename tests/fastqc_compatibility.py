"""Run FastQC and MultiQC on the CLI artifacts from a small FASTQ fixture."""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile


def run(*command, **kwargs):
    subprocess.run(command, check=True, **kwargs)


def sections(archive_path, sample):
    with ZipFile(archive_path) as archive:
        text = archive.read(f"{sample}_fastqc/fastqc_data.txt").decode("utf-8")
    data = {}
    for module in text.split(">>END_MODULE"):
        lines = module.splitlines()
        for index, line in enumerate(lines):
            if line.startswith(">>"):
                name, status = line[2:].split("\t")
                data[name] = (status, [row.split("\t") for row in lines[index + 1:]])
                break
    return data


def multiqc(root, archives):
    run(
        sys.executable, "-m", "multiqc", "-m", "fastqc", "--require-logs",
        "--no-report", "--outdir", str(root), *map(str, archives),
    )
    source = root / "multiqc_data/multiqc_fastqc.txt"
    with source.open(encoding="utf-8", newline="") as table:
        rows = list(csv.DictReader(table, delimiter="\t"))
    assert len({row["Sample"] for row in rows}) == len(rows)
    return {row["Sample"]: row for row in rows}


def verify(binary, fastqc, root):
    sample = root / "oracle.fastq"
    sequence = b"GCNN" + b"ACGT" * 24
    quality = b"I" * 9 + b"5H" + b"I" * 89
    records = [
        f"@MK:1:FLOW:1:1001:{i}:1\n".encode() + sequence + b"\n+\n" + quality + b"\n"
        for i in range(120)
    ]
    records.append(b"@MK:1:FLOW:1:1001:120:1\n" + sequence + b"A\n+\n" + quality + b"!\n")
    sample.write_bytes(b"".join(records))
    reference = root / "reference"
    reference.mkdir()
    run(
        "java", "-Djava.awt.headless=true", "-Dfastqc.quiet=true",
        f"-Dfastqc.output_dir={reference}", "-cp", f"{fastqc}{os.pathsep}{fastqc / '*'}",
        "uk.ac.babraham.FastQC.FastQCApplication", str(sample),
    )
    candidate = root / "candidate"
    options = ["--dup-window", "10", "--noextract", "--json", "--quiet",
               "-t", "2", "--outdir", str(candidate)]
    run(str(binary), str(sample), *options)
    for name in ["stdin_a", "stdin_b"]:
        with sample.open("rb") as stream:
            run(str(binary), "-", "--sample-name", name, *options, stdin=stream)

    reference_zip = reference / "oracle_fastqc.zip"
    candidate_zip = candidate / "oracle_fastqc.zip"
    lhs = sections(reference_zip, "oracle")
    rhs = sections(candidate_zip, "oracle")
    reference_basic = dict(lhs["Basic Statistics"][1][1:])
    candidate_basic = dict(rhs["Basic Statistics"][1][1:])
    for key in ["Filename", "Total Sequences", "Sequence length", "%GC"]:
        assert reference_basic[key] == candidate_basic[key], (key, reference_basic, candidate_basic)

    reference_status, reference_rows = lhs["Per base sequence quality"]
    candidate_status, candidate_rows = rhs["Per base sequence quality"]
    assert reference_status == candidate_status
    assert len(reference_rows) == len(candidate_rows)
    for expected, actual in zip(reference_rows[1:], candidate_rows[1:], strict=True):
        assert expected[0] == actual[0], (expected, actual)
        for a, b in zip(map(float, expected[1:]), map(float, actual[1:]), strict=True):
            assert (math.isnan(a) and math.isnan(b)) or math.isclose(a, b, abs_tol=0.005), (a, b)

    archives = sorted(candidate.glob("*_fastqc.zip"))
    assert len(archives) == 3
    imported = multiqc(root / "multiqc-kira", archives)
    assert set(imported) == {"oracle", "stdin_a", "stdin_b"}, imported.keys()
    for row in imported.values():
        assert float(row["Total Sequences"]) == len(records), row
        assert float(row["%GC"]) == float(reference_basic["%GC"]), row
        assert float(row["total_deduplicated_percentage"]) == 10.0, row
    assert set(multiqc(root / "multiqc-fastqc", [reference_zip])) == {"oracle"}
    with ZipFile(candidate_zip) as archive:
        report = json.loads(archive.read("oracle_fastqc/kira_qc.json"))
        assert report["duplication"]["reads_analysed"] == 10
        assert report["duplication"]["whole_file"] is False
        assert (candidate / "oracle_fastqc.html").read_bytes() == archive.read(
            "oracle_fastqc/fastqc_report.html"
        )
    print("FastQC numeric checks and MultiQC import passed for all three candidate samples")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("fastqc_directory", type=Path)
    args = parser.parse_args()
    with TemporaryDirectory() as temp:
        verify(args.binary.resolve(), args.fastqc_directory.resolve(), Path(temp))
