#!/usr/bin/env python3
"""Disjoint-cohort rustroo/Panaroo curve for one MACS allocation.

The partial cohorts are different assemblies. The last point reruns every
assembly. Rustroo and Panaroo run one after the other, so each keeps the
whole CPU count and memory of the allocation. Capture is a reciprocal BLAST
hit at >=70% identity and >=90% coverage of both reference sequences, the
same acceptance rule as the local chromosome benchmark.

The compute-node python3 can be older than 3.7, so this file stays on that
syntax and does not import the newer benchmark module.
"""

import argparse
import csv
import fcntl
import json
import math
import os
from pathlib import Path
import random
import resource
import shutil
import statistics
import subprocess
import sys
import time


CURVE_STEPS = 20
CURVE_START = 1
CURVE_STOP = 1500
MIN_IDENTITY = 70.0
MIN_COVERAGE = 90.0
EVALUE = "1e-10"
GFF_SUFFIXES = (".gff", ".gff3", ".gff.gz", ".gff3.gz")
BLAST_FORMAT = (
    "6 qseqid sseqid pident length nident qlen slen "
    "qstart qend sstart send evalue bitscore"
)


def save_json(path, value):
    # Each writer uses its own temporary file. Twenty array tasks share this
    # directory, and one shared name makes replace lose the source.
    temporary = path.with_name("{0}.{1}.tmp".format(path.name, os.getpid()))
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def load_json(path):
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError("Expected a JSON object: {0}".format(path))
    return value


def write_tsv(path, rows):
    if not rows:
        raise ValueError("No rows to write: {0}".format(path))
    fields = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    temporary = path.with_name("{0}.{1}.tmp".format(path.name, os.getpid()))
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def gff_paths(directory):
    files = [
        path for path in directory.iterdir()
        if path.is_file() and path.name.lower().endswith(GFF_SUFFIXES)
    ]
    if not files:
        raise ValueError("No GFF files in {0}".format(directory))
    return sorted(files, key=lambda path: path.name)


def curve_sizes(steps=CURVE_STEPS, start=CURVE_START, stop=CURVE_STOP):
    """`steps` sizes from `start` to `stop`. Growth is exponential, so most steps are small."""
    if steps < 2 or start < 1 or stop < start:
        raise ValueError("The curve needs at least two increasing sizes")
    log_start = math.log(start)
    log_stop = math.log(stop)
    sizes = []
    for index in range(steps):
        grown = math.exp(log_start + (log_stop - log_start) * index / float(steps - 1))
        size = int(round(grown))
        if sizes and size <= sizes[-1]:
            size = sizes[-1] + 1
        sizes.append(size)
    if sizes[-1] != stop:
        sizes[-1] = stop
    if sizes[-1] <= sizes[-2]:
        raise ValueError("The curve ends are too close for {0} steps".format(steps))
    return tuple(sizes)


def split_cohorts(names, sizes=None, seed=1729):
    """Shuffle once. Small cohorts take unused genomes. A size that no longer fits reuses them."""
    if sizes is None:
        sizes = curve_sizes()
    ordered = tuple(sorted(names))
    limit = len(ordered)
    if limit < 1:
        raise ValueError("No assemblies to split")
    resolved = []
    seen_size = set()
    for size in sizes:
        if size < 1:
            raise ValueError("Cohort sizes must be positive")
        capped = size if size <= limit else limit
        if capped in seen_size:
            continue
        seen_size.add(capped)
        resolved.append(capped)
    draw = list(ordered)
    random.Random(seed).shuffle(draw)
    cohorts = []
    cursor = 0
    for size in resolved:
        if cursor + size <= limit:
            group = tuple(draw[cursor:cursor + size])
            cursor += size
            disjoint = True
        else:
            group = tuple(draw[:size])
            disjoint = False
        if len(group) != size:
            raise RuntimeError("Cohort split failed")
        cohorts.append({"n": size, "disjoint": disjoint, "files": group})
    return cohorts


def choose_cohorts(cohorts, only=None, task=None):
    """`only` is a sample count or ``full``. ``task`` is one curve step."""
    if task is not None:
        if task < 0 or task >= len(cohorts):
            raise SystemExit("array task {0} has no cohort".format(task))
        return [cohorts[task]]
    if only is None:
        return cohorts
    if only == "full":
        chosen = [max(cohorts, key=lambda cohort: cohort["n"])]
    else:
        try:
            size = int(only)
        except (TypeError, ValueError):
            raise SystemExit("--n must be a sample count or full")
        chosen = [cohort for cohort in cohorts if cohort["n"] == size]
    if len(chosen) != 1:
        raise SystemExit("No single cohort matches --n {0}".format(only))
    return chosen


def capture_fields(metrics):
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


def fasta_sequences(path):
    sequences = {}
    identifier = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()
                if not header:
                    raise ValueError("Empty FASTA header: {0}".format(path))
                identifier = header[0]
                if identifier in sequences:
                    raise ValueError("Duplicate FASTA identifier {0}: {1}".format(identifier, path))
                sequences[identifier] = ""
            elif identifier is None:
                raise ValueError("Sequence before FASTA header: {0}".format(path))
            else:
                sequences[identifier] += line.upper()
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValueError("Empty FASTA or empty sequence: {0}".format(path))
    return sequences


def best_hits(path):
    best = {}
    with path.open() as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != 13:
                raise ValueError("Expected 13 BLAST columns: {0}:{1}".format(path, line_number))
            qlen, slen = int(row[5]), int(row[6])
            qspan = abs(int(row[8]) - int(row[7])) + 1
            sspan = abs(int(row[10]) - int(row[9])) + 1
            if (
                qlen <= 0 or slen <= 0 or qspan > qlen or sspan > slen
                or not 1 <= int(row[7]) <= qlen or not 1 <= int(row[8]) <= qlen
                or not 1 <= int(row[9]) <= slen or not 1 <= int(row[10]) <= slen
                or not 0 <= int(row[4]) <= int(row[3]) or int(row[3]) <= 0
            ):
                raise ValueError("Invalid BLAST lengths: {0}:{1}".format(path, line_number))
            score = float(row[12])
            if not math.isfinite(score) or not 0 <= float(row[2]) <= 100:
                raise ValueError("Invalid BLAST score: {0}:{1}".format(path, line_number))
            query, subject = row[:2]
            previous = best.get(query)
            if previous is None or score > previous[0]:
                best[query] = (score, {subject: row})
            elif score == previous[0]:
                old = previous[1].get(subject)

                def rank(hit):
                    coverage = min(
                        (abs(int(hit[8]) - int(hit[7])) + 1) / float(hit[5]),
                        (abs(int(hit[10]) - int(hit[9])) + 1) / float(hit[6]),
                    )
                    return coverage, float(hit[2]), tuple(hit)

                if old is None or rank(row) > rank(old):
                    previous[1][subject] = row
    ambiguous = set(query for query, (_, hits) in best.items() if len(hits) > 1)
    unique = dict(
        (query, next(iter(hits.values())))
        for query, (_, hits) in best.items() if len(hits) == 1
    )
    return unique, ambiguous, set(best)


def compare_references(forward_path, reverse_path, panaroo, rustroo, min_identity, min_coverage):
    forward, ambiguous_pan, hits_pan = best_hits(forward_path)
    reverse, ambiguous_rust, hits_rust = best_hits(reverse_path)
    if not hits_pan.issubset(set(panaroo)) or not hits_rust.issubset(set(rustroo)):
        raise ValueError("BLAST query identifiers do not match the references")
    for hits, targets, queries in ((forward, rustroo, panaroo), (reverse, panaroo, rustroo)):
        for query, row in hits.items():
            if row[1] not in targets or int(row[5]) != len(queries[query]) or int(row[6]) != len(targets[row[1]]):
                raise ValueError("BLAST targets or lengths do not match the references")
    pairs = []
    reciprocal = 0
    for query, row in sorted(forward.items()):
        subject = row[1]
        back = reverse.get(subject)
        if back is None or back[1] != query:
            continue
        reciprocal += 1
        identity = float(row[2])
        qcoverage = 100 * (abs(int(row[8]) - int(row[7])) + 1) / float(row[5])
        scoverage = 100 * (abs(int(row[10]) - int(row[9])) + 1) / float(row[6])
        if identity < min_identity or min(qcoverage, scoverage) < min_coverage:
            continue
        nident, columns = int(row[4]), int(row[3])
        exact = nident == columns == int(row[5]) == int(row[6]) and min(qcoverage, scoverage) == 100
        pairs.append({
            "panaroo_gene": query, "rustroo_gene": subject, "identity_pct": identity,
            "panaroo_coverage_pct": qcoverage, "rustroo_coverage_pct": scoverage,
            "identical_bases": nident, "alignment_columns": columns,
            "full_length_exact": int(exact),
        })
    metrics = {
        "panaroo_references": len(panaroo), "rustroo_references": len(rustroo),
        "panaroo_no_hit": len(panaroo) - len(hits_pan),
        "rustroo_no_hit": len(rustroo) - len(hits_rust),
        "panaroo_ambiguous": len(ambiguous_pan), "rustroo_ambiguous": len(ambiguous_rust),
        "panaroo_unique_nonreciprocal": len(forward) - reciprocal,
        "rustroo_unique_nonreciprocal": len(reverse) - reciprocal,
        "reciprocal_pairs": reciprocal, "reciprocal_filtered_out": reciprocal - len(pairs),
        "accepted_pairs": len(pairs), "full_length_exact": sum(pair["full_length_exact"] for pair in pairs),
        "pairs_at_least_99pct_identity": sum(pair["identity_pct"] >= 99 for pair in pairs),
        "panaroo_accepted_fraction": len(pairs) / float(len(panaroo)),
        "rustroo_accepted_fraction": len(pairs) / float(len(rustroo)),
        "mean_identity_pct": statistics.mean([pair["identity_pct"] for pair in pairs]) if pairs else None,
        "median_identity_pct": statistics.median([pair["identity_pct"] for pair in pairs]) if pairs else None,
        "weighted_identity_pct": (
            100.0 * sum(pair["identical_bases"] for pair in pairs) / sum(pair["alignment_columns"] for pair in pairs)
        ) if pairs else None,
    }
    return metrics, pairs


def _minutes(seconds):
    return seconds / 60.0


def run_timed(argv, log_path, env):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    with log_path.open("w") as handle:
        completed = subprocess.run(argv, stdout=handle, stderr=subprocess.STDOUT, env=env)
    wall = time.perf_counter() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    return completed.returncode, wall, max(cpu, 0.0)


def reported_rustroo_minutes(log_path):
    found = {}
    if log_path.is_file():
        with log_path.open() as handle:
            for line in handle:
                key, separator, value = line.rstrip("\n").partition("\t")
                if separator and key in ("cpu_min", "actual_min"):
                    try:
                        found[key] = float(value)
                    except ValueError:
                        continue
    return {
        "rustroo_reported_cpu_min": found.get("cpu_min"),
        "rustroo_reported_actual_min": found.get("actual_min"),
    }


def rustroo_argv(binary, gffs, output, threads, quality="standard"):
    return [str(binary), "-i"] + [str(path) for path in gffs] + [
        "-o", str(output), "--threads", str(threads), "--quality", quality,
    ]


def panaroo_argv(sif, gffs, cohort_dir, threads, input_root):
    # Bind the cohort directory, which already exists. Panaroo creates
    # panaroo/ inside it; binding that subdirectory first would require it
    # to exist before Panaroo is allowed to create it.
    inside = ["/inputs/{0}".format(path.name) for path in gffs]
    # One argument is read as a list of paths, not as a GFF. A real GFF then
    # fails with an empty extension. Hand Panaroo a one-line list instead.
    if len(inside) == 1:
        (cohort_dir / "panaroo_inputs.txt").write_text("{0}\n".format(inside[0]))
        inside = ["/cohort/panaroo_inputs.txt"]
    return [
        "apptainer", "exec", "--cleanenv", "--no-home",
        "--bind", "{0}:/inputs:ro".format(input_root),
        "--bind", "{0}:/cohort".format(cohort_dir),
        str(sif),
        "env", "OMP_NUM_THREADS={0}".format(threads), "PYTHONHASHSEED=0",
        "panaroo", "-i",
    ] + inside + [
        "-o", "/cohort/panaroo",
        "--clean-mode", "strict", "--remove-invalid-genes",
        "-t", str(threads),
    ]


def blast_argv(sif, work, threads, rust_count, pan_count):
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
        "--bind", "{0}:/work".format(work), str(sif),
        "bash", "-euc", script, "bash",
        str(threads), EVALUE, str(rust_count), str(pan_count), BLAST_FORMAT,
    ]


def cohort_folder(n):
    return "n{0:04d}".format(n)


def report_failure(row, result_path, message, log_path):
    row["error"] = message
    save_json(result_path, row)
    print(message, flush=True)
    if log_path.is_file():
        with log_path.open(errors="replace") as handle:
            lines = handle.read().splitlines()
        if lines:
            print("\n".join(lines[-20:]), flush=True)
    return row


def run_cohort(root, cohort, input_root, binary, sif, threads, env, quality):
    folder = root / cohort_folder(cohort["n"])
    result_path = folder / "result.json"
    if result_path.is_file() and load_json(result_path).get("status") == "ok":
        print("skip {0}: already complete".format(folder.name), flush=True)
        return load_json(result_path)
    if folder.exists():
        folder.rename(folder.with_name("{0}.failed_{1}".format(
            folder.name, time.strftime("%Y%m%d_%H%M%S"),
        )))
    folder.mkdir(parents=True)
    names = cohort["files"]
    gffs = [input_root / name for name in names]
    missing = [str(path) for path in gffs if not path.is_file()]
    if missing:
        raise ValueError("Cohort file missing: {0}".format(missing[0]))
    (folder / "inputs.txt").write_text("".join("{0}\n".format(name) for name in names))
    row = {
        "n": cohort["n"],
        "disjoint": cohort["disjoint"],
        "gff_count": len(names),
        "threads": threads,
        "quality": quality,
        "status": "failed",
        "error": "",
    }
    save_json(result_path, row)
    print("cohort n={0} disjoint={1}".format(cohort["n"], cohort["disjoint"]), flush=True)

    # One directory per cohort, next to the references, so a download of the
    # run folder includes every tool log.
    logs = folder / "logs"
    logs.mkdir()
    print("logs {0}".format(logs), flush=True)
    rust_dir = folder / "rustroo"
    code, wall, cpu = run_timed(
        rustroo_argv(binary, gffs, rust_dir, threads, quality), logs / "rustroo.stdout", env,
    )
    detail = rust_dir / "rustroo.log"
    if detail.is_file():
        shutil.copy2(detail, logs / "rustroo.log")
    row.update({
        "rustroo_exit": code,
        "rustroo_wall_min": _minutes(wall),
        "rustroo_cpu_min": _minutes(cpu),
    })
    row.update(reported_rustroo_minutes(logs / "rustroo.stdout"))
    save_json(result_path, row)
    if code != 0:
        return report_failure(
            row, result_path, "rustroo exited {0}".format(code), logs / "rustroo.stdout",
        )

    pan_dir = folder / "panaroo"
    code, wall, cpu = run_timed(
        panaroo_argv(sif, gffs, folder, threads, input_root), logs / "panaroo.stdout", env,
    )
    row.update({
        "panaroo_exit": code,
        "panaroo_wall_min": _minutes(wall),
        "panaroo_cpu_min": _minutes(cpu),
    })
    save_json(result_path, row)
    if code != 0:
        return report_failure(
            row, result_path, "panaroo exited {0}".format(code), logs / "panaroo.stdout",
        )

    rust_fa = rust_dir / "genes.fa"
    pan_fa = pan_dir / "pan_genome_reference.fa"
    if not rust_fa.is_file() or not pan_fa.is_file():
        return report_failure(
            row, result_path, "A reference FASTA is missing after both tools exited 0",
            logs / "panaroo.stdout",
        )
    rust_seqs = fasta_sequences(rust_fa)
    pan_seqs = fasta_sequences(pan_fa)
    row.update({"rustroo_genes": len(rust_seqs), "panaroo_genes": len(pan_seqs)})
    work = folder / "compare"
    work.mkdir()
    shutil.copy2(pan_fa, work / "panaroo.fa")
    shutil.copy2(rust_fa, work / "rustroo.fa")
    code, wall, cpu = run_timed(
        blast_argv(sif, work, threads, len(rust_seqs), len(pan_seqs)),
        logs / "blast.stdout", env,
    )
    row.update({
        "blast_exit": code,
        "blast_wall_min": _minutes(wall),
        "blast_cpu_min": _minutes(cpu),
    })
    if code != 0:
        return report_failure(
            row, result_path, "blastn exited {0}".format(code), logs / "blast.stdout",
        )
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
    row.update(capture_fields(metrics))
    row["status"] = "ok"
    row["error"] = ""
    save_json(result_path, row)
    identity = row["mean_identity_pct"]
    identity_text = "{0:.2f}".format(identity) if identity is not None else "none"
    print(
        "n={0} rustroo {1:.2f} min, panaroo {2:.2f} min, panaroo captured {3:.1f}%, rustroo captured {4:.1f}%, mean identity {5}".format(
            row["n"], row["rustroo_wall_min"], row["panaroo_wall_min"],
            row["panaroo_captured_pct"], row["rustroo_captured_pct"], identity_text,
        ),
        flush=True,
    )
    return row


def _series(rows, key):
    points = []
    for row in rows:
        value = row.get(key)
        if row.get("status") == "ok" and isinstance(value, (int, float)):
            points.append((int(row["n"]), float(value)))
    return points


def plot(path, title, ylabel, series):
    present = dict((name, points) for name, points in series.items() if points)
    if not present:
        return
    sizes = sorted(set(n for points in present.values() for n, _ in points))
    maximum = max(value for points in present.values() for _, value in points) * 1.1 or 1.0
    left, top, width, height = 90, 55, 660, 300
    span = (max(sizes) - min(sizes)) or 1

    def x_of(n):
        return left + (n - min(sizes)) / float(span) * width

    def y_of(value):
        return top + height * (1 - value / maximum)

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="430" viewBox="0 0 900 430">',
        '<rect width="900" height="430" fill="white"/>',
        '<g font-family="sans-serif" font-size="13">',
        '<text x="90" y="28" font-size="18">{0}</text>'.format(_xml(title)),
        '<text x="420" y="408" text-anchor="middle">{0}</text>'.format(_xml(ylabel)),
    ]
    for tick in range(6):
        value = maximum * tick / 5.0
        yy = y_of(value)
        lines.append(
            '<line x1="{0}" y1="{1:.1f}" x2="{2}" y2="{1:.1f}" stroke="#ddd"/>'.format(left, yy, left + width)
        )
        lines.append('<text x="82" y="{0:.1f}" text-anchor="end">{1:.3g}</text>'.format(yy + 4, value))
    labeled_at = None
    for n in sizes:
        x = x_of(n)
        if labeled_at is not None and x - labeled_at < 46 and n != sizes[-1]:
            continue
        if labeled_at is not None and x - labeled_at < 46:
            lines.pop()
        lines.append('<text x="{0:.1f}" y="377" text-anchor="middle">{1}</text>'.format(x, n))
        labeled_at = x
    lines.append('<text x="420" y="396" text-anchor="middle">Assemblies</text>')
    for index, name in enumerate(present):
        points = present[name]
        color = ("#0072B2", "#D55E00", "#009E73")[index % 3]
        ordered = sorted(points)
        coords = " ".join("{0:.1f},{1:.1f}".format(x_of(n), y_of(value)) for n, value in ordered)
        lines.append('<polyline points="{0}" fill="none" stroke="{1}" stroke-width="2"/>'.format(coords, color))
        for n, value in ordered:
            lines.append('<circle cx="{0:.1f}" cy="{1:.1f}" r="4" fill="{2}"/>'.format(x_of(n), y_of(value), color))
        lines.append('<text x="770" y="{0}" fill="{1}">{2}</text>'.format(80 + index * 22, color, _xml(name)))
    lines.append("</g></svg>\n")
    temporary = path.with_name("{0}.{1}.tmp".format(path.name, os.getpid()))
    temporary.write_text("\n".join(lines))
    temporary.replace(path)


def _xml(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def publish(root):
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".publish.lock").open("a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        return _publish(root)


def _publish(root):
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


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rustroo", type=Path, required=True)
    parser.add_argument("--panaroo-sif", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--quality", default="standard", choices=("standard", "high", "strict"))
    parser.add_argument("--n", default=None, help="one sample count, or full, on this node")
    parser.add_argument("--task", type=int, default=None, help="one curve step, counting from 0")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    if args.threads < 1:
        raise SystemExit("threads must be positive")
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None and int(slurm_cpus) != args.threads:
        raise SystemExit(
            "--threads {0} does not match SLURM_CPUS_PER_TASK={1}".format(args.threads, slurm_cpus)
        )
    if not args.rustroo.is_file() or not args.panaroo_sif.is_file():
        raise SystemExit("rustroo binary or Panaroo image is missing")
    names = [path.name for path in gff_paths(args.input)]
    cohorts = choose_cohorts(
        split_cohorts(names, curve_sizes(), args.seed), args.n, args.task,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    if args.task is not None:
        manifest = "run.task{0}.json".format(args.task)
    elif args.n is None:
        manifest = "run.json"
    else:
        manifest = "run.{0}.json".format(args.n)
    save_json(args.output / manifest, {
        "seed": args.seed,
        "threads": args.threads,
        "quality": args.quality,
        "slurm_cpus": slurm_cpus,
        "slurm_mem_mb": os.environ.get("SLURM_MEM_PER_NODE"),
        "partial_sizes": [cohort["n"] for cohort in cohorts],
        "full_n": len(names),
        "identity_pct": MIN_IDENTITY,
        "coverage_pct": MIN_COVERAGE,
        "n": args.n,
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
                args.output, cohort, args.input, args.rustroo, args.panaroo_sif,
                args.threads, env, args.quality,
            )
            failed += row.get("status") != "ok"
            publish(args.output)
    finally:
        publish(args.output)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
