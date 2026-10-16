#!/usr/bin/env python3
"""Disjoint-cohort rustroo/Panaroo curve for one MACS allocation.

The partial cohorts are different assemblies. The last point reruns every
assembly. Rustroo and Panaroo run one after the other, so each keeps the
whole CPU count and memory of the allocation. Capture is a reciprocal BLAST
hit at >=70% identity and >=90% coverage of both reference sequences, the
same acceptance rule as the local chromosome benchmark.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import random
import resource
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable, Optional

from benchmark import (
    BLAST_FORMAT,
    compare_references,
    fasta_sequences,
    load_json,
    save_json,
    write_tsv,
)


PARTIAL_SIZES = (50, 150, 300, 800)
MIN_IDENTITY = 70.0
MIN_COVERAGE = 90.0
EVALUE = "1e-10"
GFF_SUFFIXES = (".gff", ".gff3", ".gff.gz", ".gff3.gz")


def gff_paths(directory: Path) -> list[Path]:
    files = [
        path for path in directory.iterdir()
        if path.is_file() and path.name.lower().endswith(GFF_SUFFIXES)
    ]
    if not files:
        raise ValueError(f"No GFF files in {directory}")
    return sorted(files, key=lambda path: path.name)


def split_cohorts(
    names: Iterable[str],
    sizes: tuple[int, ...] = PARTIAL_SIZES,
    seed: int = 1729,
) -> list[dict[str, Any]]:
    """Shuffle once, then cut non-overlapping cohorts. The full list is last."""
    ordered = tuple(sorted(names))
    if any(size < 1 for size in sizes):
        raise ValueError("Cohort sizes must be positive")
    needed = sum(sizes)
    if needed > len(ordered):
        raise ValueError(
            f"Disjoint cohorts need {needed} assemblies; the input has {len(ordered)}"
        )
    draw = list(ordered)
    random.Random(seed).shuffle(draw)
    cohorts: list[dict[str, Any]] = []
    seen: set[str] = set()
    cursor = 0
    for size in sizes:
        group = tuple(draw[cursor:cursor + size])
        cursor += size
        if len(group) != size or seen.intersection(group):
            raise RuntimeError("Disjoint split failed")
        seen.update(group)
        cohorts.append({"n": size, "disjoint": True, "files": group})
    cohorts.append({"n": len(ordered), "disjoint": False, "files": ordered})
    return cohorts


def capture_fields(metrics: dict[str, Any]) -> dict[str, Any]:
    """Percent of each tool's references with an accepted match in the other."""
    pan_total = metrics["panaroo_references"]
    rust_total = metrics["rustroo_references"]
    matched = metrics["accepted_pairs"]
    if pan_total < 1 or rust_total < 1:
        raise ValueError("Capture requires at least one reference from each tool")
    return {
        "panaroo_captured_pct": 100.0 * matched / pan_total,
        "rustroo_captured_pct": 100.0 * matched / rust_total,
        "mean_identity_pct": metrics["mean_identity_pct"],
        "weighted_identity_pct": metrics["weighted_identity_pct"],
        "full_length_exact": metrics["full_length_exact"],
        "accepted_pairs": matched,
    }


def _minutes(seconds: float) -> float:
    return seconds / 60.0


def run_timed(argv: list[str], log_path: Path, env: dict[str, str]) -> tuple[int, float, float]:
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    with log_path.open("w") as handle:
        completed = subprocess.run(
            argv, stdout=handle, stderr=subprocess.STDOUT, env=env, check=False,
        )
    wall = time.perf_counter() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    return completed.returncode, wall, max(cpu, 0.0)


def reported_rustroo_minutes(log_path: Path) -> dict[str, Optional[float]]:
    found: dict[str, float] = {}
    if log_path.is_file():
        for line in log_path.read_text(errors="replace").splitlines():
            key, separator, value = line.partition("\t")
            if separator and key in {"cpu_min", "actual_min"}:
                try:
                    found[key] = float(value)
                except ValueError:
                    continue
    return {
        "rustroo_reported_cpu_min": found.get("cpu_min"),
        "rustroo_reported_actual_min": found.get("actual_min"),
    }


def rustroo_argv(binary: Path, gffs: Iterable[Path], output: Path, threads: int) -> list[str]:
    return [
        str(binary), "-i", *[str(path) for path in gffs],
        "-o", str(output), "--threads", str(threads),
    ]


def panaroo_argv(
    sif: Path, gffs: Iterable[Path], cohort_dir: Path, threads: int, input_root: Path,
) -> list[str]:
    # Bind the cohort directory, which already exists. Panaroo creates
    # `panaroo/` inside it; binding that subdirectory first would require it
    # to exist before Panaroo is allowed to create it.
    inside = [f"/inputs/{path.name}" for path in gffs]
    return [
        "apptainer", "exec", "--cleanenv", "--no-home",
        "--bind", f"{input_root}:/inputs:ro",
        "--bind", f"{cohort_dir}:/cohort",
        str(sif),
        "env", f"OMP_NUM_THREADS={threads}", "PYTHONHASHSEED=0",
        "panaroo", "-i", *inside, "-o", "/cohort/panaroo",
        "--clean-mode", "strict", "--remove-invalid-genes", "-t", str(threads),
    ]


def blast_argv(sif: Path, work: Path, threads: int, rust_count: int, pan_count: int) -> list[str]:
    script = """
set -euo pipefail
makeblastdb -in /work/panaroo.fa -dbtype nucl -out /work/panaroo_db
makeblastdb -in /work/rustroo.fa -dbtype nucl -out /work/rustroo_db
blastn -task blastn -query /work/panaroo.fa -db /work/rustroo_db -num_threads "$1" \
  -evalue "$2" -dust no -max_target_seqs "$3" -outfmt "$5" -out /work/forward.tsv
blastn -task blastn -query /work/rustroo.fa -db /work/panaroo_db -num_threads "$1" \
  -evalue "$2" -dust no -max_target_seqs "$4" -outfmt "$5" -out /work/reverse.tsv
"""
    return [
        "apptainer", "exec", "--cleanenv", "--no-home",
        "--bind", f"{work}:/work", str(sif),
        "bash", "-euc", script, "bash",
        str(threads), EVALUE, str(rust_count), str(pan_count), BLAST_FORMAT,
    ]


def cohort_folder(n: int) -> str:
    return f"n{n:04d}"


def run_cohort(
    root: Path,
    cohort: dict[str, Any],
    input_root: Path,
    binary: Path,
    sif: Path,
    threads: int,
    env: dict[str, str],
) -> dict[str, Any]:
    folder = root / cohort_folder(cohort["n"])
    result_path = folder / "result.json"
    if result_path.is_file() and load_json(result_path).get("status") == "ok":
        print(f"skip {folder.name}: already complete", flush=True)
        return load_json(result_path)
    if folder.exists():
        folder.rename(folder.with_name(f"{folder.name}.failed_{time.strftime('%Y%m%d_%H%M%S')}"))
    folder.mkdir(parents=True)
    names = cohort["files"]
    gffs = [input_root / name for name in names]
    missing = [str(path) for path in gffs if not path.is_file()]
    if missing:
        raise ValueError(f"Cohort file missing: {missing[0]}")
    (folder / "inputs.txt").write_text("".join(f"{name}\n" for name in names))
    row: dict[str, Any] = {
        "n": cohort["n"],
        "disjoint": cohort["disjoint"],
        "gff_count": len(names),
        "threads": threads,
        "status": "failed",
        "error": "",
    }
    save_json(result_path, row)
    print(f"cohort n={cohort['n']} disjoint={cohort['disjoint']}", flush=True)

    rust_dir = folder / "rustroo"
    code, wall, cpu = run_timed(
        rustroo_argv(binary, gffs, rust_dir, threads), folder / "rustroo.stdout", env,
    )
    row.update(
        rustroo_exit=code,
        rustroo_wall_min=_minutes(wall),
        rustroo_cpu_min=_minutes(cpu),
        **reported_rustroo_minutes(folder / "rustroo.stdout"),
    )
    save_json(result_path, row)
    if code != 0:
        row["error"] = f"rustroo exited {code}"
        save_json(result_path, row)
        return row

    pan_dir = folder / "panaroo"
    code, wall, cpu = run_timed(
        panaroo_argv(sif, gffs, folder, threads, input_root), folder / "panaroo.stdout", env,
    )
    row.update(panaroo_exit=code, panaroo_wall_min=_minutes(wall), panaroo_cpu_min=_minutes(cpu))
    save_json(result_path, row)
    if code != 0:
        row["error"] = f"panaroo exited {code}"
        save_json(result_path, row)
        return row

    rust_fa = rust_dir / "genes.fa"
    pan_fa = pan_dir / "pan_genome_reference.fa"
    if not rust_fa.is_file() or not pan_fa.is_file():
        row["error"] = "A reference FASTA is missing after both tools exited 0"
        save_json(result_path, row)
        return row
    rust_seqs = fasta_sequences(rust_fa)
    pan_seqs = fasta_sequences(pan_fa)
    row.update(rustroo_genes=len(rust_seqs), panaroo_genes=len(pan_seqs))
    work = folder / "compare"
    work.mkdir()
    shutil.copy2(pan_fa, work / "panaroo.fa")
    shutil.copy2(rust_fa, work / "rustroo.fa")
    code, wall, cpu = run_timed(
        blast_argv(sif, work, threads, len(rust_seqs), len(pan_seqs)),
        folder / "blast.stdout", env,
    )
    row.update(blast_exit=code, blast_wall_min=_minutes(wall), blast_cpu_min=_minutes(cpu))
    if code != 0:
        row["error"] = f"blastn exited {code}"
        save_json(result_path, row)
        return row
    metrics, pairs = compare_references(
        work / "forward.tsv", work / "reverse.tsv", pan_seqs, rust_seqs,
        MIN_IDENTITY, MIN_COVERAGE,
    )
    if pairs:
        write_tsv(work / "shared_genes.tsv", pairs)
    else:
        (work / "shared_genes.tsv").write_text(
            "panaroo_gene\trustroo_gene\tidentity_pct\tpanaroo_coverage_pct\t"
            "rustroo_coverage_pct\tidentical_bases\talignment_columns\tfull_length_exact\n"
        )
    row.update(capture_fields(metrics), status="ok", error="")
    save_json(result_path, row)
    print(
        f"n={row['n']} rustroo {row['rustroo_wall_min']:.2f} min, "
        f"panaroo {row['panaroo_wall_min']:.2f} min, "
        f"panaroo captured {row['panaroo_captured_pct']:.1f}%, "
        f"rustroo captured {row['rustroo_captured_pct']:.1f}%, "
        f"mean identity {row['mean_identity_pct']}",
        flush=True,
    )
    return row


def _series(rows: list[dict[str, Any]], key: str) -> list[tuple[int, float]]:
    points = []
    for row in rows:
        value = row.get(key)
        if row.get("status") == "ok" and isinstance(value, (int, float)):
            points.append((int(row["n"]), float(value)))
    return points


def plot(path: Path, title: str, ylabel: str, series: dict[str, list[tuple[int, float]]]) -> None:
    present = {name: points for name, points in series.items() if points}
    if not present:
        return
    sizes = sorted({n for points in present.values() for n, _ in points})
    maximum = max(value for points in present.values() for _, value in points) * 1.1 or 1.0
    left, top, width, height = 90, 55, 660, 300
    span = (max(sizes) - min(sizes)) or 1

    def x(n: int) -> float:
        return left + (n - min(sizes)) / span * width

    def y(value: float) -> float:
        return top + height * (1 - value / maximum)

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="430" viewBox="0 0 900 430">',
        '<rect width="900" height="430" fill="white"/>',
        '<g font-family="sans-serif" font-size="13">',
        f'<text x="90" y="28" font-size="18">{_xml(title)}</text>',
        f'<text x="420" y="408" text-anchor="middle">{_xml(ylabel)}</text>',
    ]
    for tick in range(6):
        value = maximum * tick / 5
        lines.append(
            f'<line x1="{left}" y1="{y(value):.1f}" x2="{left + width}" y2="{y(value):.1f}" stroke="#ddd"/>'
        )
        lines.append(f'<text x="82" y="{y(value) + 4:.1f}" text-anchor="end">{value:.3g}</text>')
    for n in sizes:
        lines.append(f'<text x="{x(n):.1f}" y="377" text-anchor="middle">{n}</text>')
    lines.append('<text x="420" y="396" text-anchor="middle">Assemblies</text>')
    for index, (name, points) in enumerate(present.items()):
        color = ("#0072B2", "#D55E00", "#009E73")[index % 3]
        ordered = sorted(points)
        coords = " ".join(f"{x(n):.1f},{y(value):.1f}" for n, value in ordered)
        lines.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2"/>')
        for n, value in ordered:
            lines.append(f'<circle cx="{x(n):.1f}" cy="{y(value):.1f}" r="4" fill="{color}"/>')
        lines.append(f'<text x="770" y="{80 + index * 22}" fill="{color}">{_xml(name)}</text>')
    lines.append("</g></svg>\n")
    path.write_text("\n".join(lines))


def _xml(text: str) -> str:
    return (
        text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )


def publish(root: Path) -> list[dict[str, Any]]:
    rows = []
    for result_path in sorted(root.glob("n[0-9][0-9][0-9][0-9]/result.json")):
        row = load_json(result_path)
        if "n" in row:
            rows.append(row)
    rows.sort(key=lambda row: (0 if row.get("disjoint") else 1, int(row["n"])))
    if not rows:
        return rows
    write_tsv(root / "curve.tsv", rows)
    plot(root / "curve_time.svg", "Runtime", "Wall minutes", {
        "rustroo": _series(rows, "rustroo_wall_min"),
        "panaroo": _series(rows, "panaroo_wall_min"),
    })
    plot(root / "curve_accuracy.svg", "Reference agreement", "Percent", {
        "panaroo genes in rustroo": _series(rows, "panaroo_captured_pct"),
        "rustroo genes in panaroo": _series(rows, "rustroo_captured_pct"),
        "mean alignment identity": _series(rows, "mean_identity_pct"),
    })
    return rows


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rustroo", type=Path, required=True)
    parser.add_argument("--panaroo-sif", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--seed", type=int, default=1729)
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.threads < 1:
        raise SystemExit("threads must be positive")
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None and int(slurm_cpus) != args.threads:
        raise SystemExit(
            f"--threads {args.threads} does not match SLURM_CPUS_PER_TASK={slurm_cpus}"
        )
    if not args.rustroo.is_file() or not args.panaroo_sif.is_file():
        raise SystemExit("rustroo binary or Panaroo image is missing")
    names = [path.name for path in gff_paths(args.input)]
    cohorts = split_cohorts(names, PARTIAL_SIZES, args.seed)
    args.output.mkdir(parents=True, exist_ok=True)
    save_json(args.output / "run.json", {
        "seed": args.seed,
        "threads": args.threads,
        "slurm_cpus": slurm_cpus,
        "slurm_mem_mb": os.environ.get("SLURM_MEM_PER_NODE"),
        "partial_sizes": list(PARTIAL_SIZES),
        "full_n": len(names),
        "identity_pct": MIN_IDENTITY,
        "coverage_pct": MIN_COVERAGE,
        "sequential": True,
        "cohorts": [
            {"n": cohort["n"], "disjoint": cohort["disjoint"], "files": list(cohort["files"])}
            for cohort in cohorts
        ],
    })
    env = os.environ.copy()
    env["RAYON_NUM_THREADS"] = str(args.threads)
    env["OMP_NUM_THREADS"] = str(args.threads)
    failed = 0
    try:
        for cohort in cohorts:
            row = run_cohort(
                args.output, cohort, args.input, args.rustroo, args.panaroo_sif, args.threads, env,
            )
            failed += row.get("status") != "ok"
            publish(args.output)
    finally:
        publish(args.output)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
