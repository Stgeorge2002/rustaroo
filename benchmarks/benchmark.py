#!/usr/bin/env python3
"""Reproducible chromosome-only rustroo/Panaroo benchmark (Python 3.9+)."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import platform
import random
import shutil
import signal
import statistics
import subprocess
import sys
import time
from typing import Any


SIZES = (3, 5, 10, 15, 20, 25, 30, 35, 40)
TOOLS = ("rustroo", "panaroo")
TIME_FORMAT = (
    "wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\n"
    "max_process_rss_kib=%M\nexit_code=%x"
)
BLAST_FORMAT = (
    "6 qseqid sseqid pident length nident qlen slen "
    "qstart qend sstart send evalue bitscore"
)
THREAD_ENV = {
    "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "LC_ALL": "C.UTF-8",
}


def save_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def capture(command: list[str]) -> str:
    try:
        return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).strip()
    except subprocess.CalledProcessError as error:
        output = (error.output or "").strip()
        raise RuntimeError(
            f"Command exited {error.returncode}: {command!r}\n"
            f"{output or 'The command produced no diagnostic output.'}"
        ) from error


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_design(
    names: list[str], seeds: list[int], sizes: list[int], repeats: int, order_seed: int,
    warmups: bool = True,
) -> tuple[dict[str, list[str]], list[dict[str, Any]]]:
    if not seeds or len(names) != len(set(names)) or len(seeds) != len(set(seeds)):
        raise ValueError("Input names and subset seeds must be unique")
    if not sizes or min(sizes) < 3 or max(sizes) > len(names):
        raise ValueError("Sample counts must be between 3 and the input count")
    if len(sizes) != len(set(sizes)) or repeats < 1:
        raise ValueError("Sample counts must be unique and repeats positive")
    subsets = {}
    blocks = []
    for series, seed in enumerate(seeds, 1):
        permutation = sorted(names)
        random.Random(seed).shuffle(permutation)
        for count in sorted(sizes):
            key = f"s{series:02d}_n{count:02d}"
            subsets[key] = sorted(permutation[:count])
            blocks.append({"subset": key, "series": series, "seed": seed, "n": count})
    rng = random.Random(order_seed)
    rng.shuffle(blocks)
    jobs = []
    for block in blocks:
        for repeat in range(0 if warmups else 1, repeats + 1):
            order = list(TOOLS)
            rng.shuffle(order)
            for tool in order:
                jobs.append({
                    **block, "repeat": repeat, "tool": tool,
                    "warmup": repeat == 0,
                    "id": f"{block['subset']}_r{repeat:02d}_{tool}",
                })
    return subsets, jobs


def docker_base(config: dict[str, Any], root: Path) -> list[str]:
    return [
        "docker", "run", "--rm", "--pull=never", "--network=none",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--cpuset-cpus", config["cpus"],
        "-v", f"{root / 'inputs'}:/inputs:ro",
    ]


def tool_command(config: dict[str, Any], root: Path, job: dict[str, Any]) -> list[str]:
    folder = root / "jobs" / job["id"]
    names = config["subsets"][job["subset"]]
    if job["tool"] == "rustroo":
        return [
            "taskset", "--cpu-list", config["cpus"],
            "/usr/bin/time", "-f", TIME_FORMAT, "-o", str(folder / "timing.txt"),
            str(root / "rustroo"), "-i",
            *[str(root / "inputs" / name) for name in names],
            "-o", str(folder / "output"), "--threads", str(config["threads"]),
            "--precluster-identity", "0.98", "--family-identity", "0.70",
            "--length-ratio", "0.95", "--min-cds-bp", "90",
        ]
    env_args = [
        argument for key, value in config["tool_environment"].items()
        for argument in ("-e", f"{key}={value}")
    ]
    return [
        *docker_base(config, root),
        "--name", f"rustroo_bench_{config['run_id']}_{job['id']}",
        "-v", f"{folder}:/job", "-w", "/job", *env_args, config["image_id"],
        "/usr/bin/time", "-f", TIME_FORMAT, "-o", "/job/timing.txt",
        "panaroo", "-i", *[f"/inputs/{name}" for name in names],
        "-o", "/job/output", "--clean-mode", config["clean_mode"],
        "--remove-invalid-genes", "-t", str(config["threads"]),
    ]


def execute(
    command: list[str], folder: Path, environment: dict[str, str],
    container_name: str | None = None,
) -> float:
    save_json(folder / "command.json", command)
    started = time.perf_counter()
    with (folder / "tool.log").open("w") as log:
        process = subprocess.Popen(
            command, stdout=log, stderr=subprocess.STDOUT,
            env={**os.environ, **environment}, start_new_session=True,
        )
        try:
            code = process.wait()
        except KeyboardInterrupt:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            if container_name:
                cleanup = subprocess.run(
                    ["docker", "rm", "-f", container_name],
                    stdout=log, stderr=subprocess.STDOUT, check=False,
                )
                if cleanup.returncode != 0:
                    print(
                        f"benchmark: container cleanup returned {cleanup.returncode}; "
                        f"check {container_name} and {folder / 'tool.log'}",
                        file=sys.stderr,
                    )
            raise
    elapsed = time.perf_counter() - started
    if code != 0:
        raise RuntimeError(f"Exit {code}; inspect {folder / 'tool.log'}")
    return elapsed


def timing_metrics(path: Path) -> dict[str, float]:
    values = {}
    for line in path.read_text().splitlines():
        key, separator, value = line.partition("=")
        if not separator or key in values:
            raise ValueError(f"Invalid GNU time output: {path}")
        number = float(value)
        if not math.isfinite(number) or number < 0:
            raise ValueError(f"Invalid timing value: {line}")
        values[key] = number
    expected = {
        "wall_seconds", "user_seconds", "system_seconds",
        "max_process_rss_kib", "exit_code",
    }
    if set(values) != expected or values["exit_code"] != 0 or values["wall_seconds"] <= 0:
        raise ValueError(f"Incomplete or unsuccessful timing: {path}")
    values["average_cpu_cores"] = (
        values["user_seconds"] + values["system_seconds"]
    ) / values["wall_seconds"]
    return values


def fasta_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    identifier = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()
                if not header:
                    raise ValueError(f"Empty FASTA header: {path}")
                identifier = header[0]
                if identifier in sequences:
                    raise ValueError(f"Duplicate FASTA identifier {identifier}: {path}")
                sequences[identifier] = ""
            elif identifier is None:
                raise ValueError(f"Sequence before FASTA header: {path}")
            else:
                sequences[identifier] += line.upper()
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValueError(f"Empty FASTA or empty sequence: {path}")
    return sequences


def output_metrics(folder: Path, tool: str, names: list[str]) -> dict[str, Any]:
    fasta = folder / ("genes.fa" if tool == "rustroo" else "pan_genome_reference.fa")
    matrix = folder / ("presence.tsv" if tool == "rustroo" else "gene_presence_absence.Rtab")
    sequences = fasta_sequences(fasta)
    groups = core = soft_core = 0
    gene_ids = set()
    with matrix.open() as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty presence matrix: {matrix}")
        samples = header.split()[1:]
        if len(samples) != len(names) or set(samples) != {Path(name).stem for name in names}:
            raise ValueError(f"Wrong samples in {matrix}")
        for line in handle:
            fields = line.split()
            if len(fields) != len(samples) + 1 or any(v not in ("0", "1") for v in fields[1:]):
                raise ValueError(f"Invalid presence row in {matrix}")
            if fields[0] in gene_ids:
                raise ValueError(f"Duplicate gene in {matrix}: {fields[0]}")
            gene_ids.add(fields[0])
            present = sum(map(int, fields[1:]))
            groups += 1
            core += 100 * present >= 99 * len(samples)
            soft_core += 95 * len(samples) <= 100 * present < 99 * len(samples)
    if not groups or not set(sequences).issubset(gene_ids):
        raise ValueError(f"References do not correspond to presence groups: {folder}")
    normalized = json.dumps(sorted(sequences.items()), separators=(",", ":")).encode()
    return {
        "groups": groups, "references": len(sequences),
        "core_groups_99pct": core, "soft_core_groups_95pct": soft_core,
        "reference_sha256": hashlib.sha256(normalized).hexdigest(),
        "presence_file_sha256": digest(matrix),
        "invalid_cds_log_messages": sum(
            "invalid gene! file - id:" in line
            for line in (folder.parent / "tool.log").read_text().splitlines()
        ),
    }


def prepare(args: argparse.Namespace) -> None:
    if sys.platform != "linux":
        raise ValueError("Run this benchmark in Linux/WSL, not Windows Python")
    source = args.input.resolve()
    binary = args.binary.resolve()
    root = args.output.resolve()
    names = sorted(path.name for path in source.glob("*.gff3") if path.is_file())
    seeds = [args.seed + index for index in range(args.subsets)]
    sizes = sorted(args.sizes)
    subsets, jobs = make_design(
        names, seeds, sizes, args.repeats, args.seed + 10000,
        warmups=not args.no_warmups,
    )
    for name in names:
        if any(character.isspace() for character in name) or ":" in name:
            raise ValueError(f"Unsupported input filename: {name}")
        with (source / name).open("rb") as handle:
            if not any(line.strip() == b"##FASTA" for line in handle):
                raise ValueError(f"Embedded ##FASTA required: {source / name}")
    with binary.open("rb") as handle:
        if handle.read(4) != b"\x7fELF":
            raise ValueError("The rustroo binary must be a Linux ELF executable")
    available = sorted(os.sched_getaffinity(0))
    threads = len(available) if args.threads is None else args.threads
    if threads < 1 or threads > len(available):
        raise ValueError(f"Threads must be in 1..{len(available)}")
    if args.subsets < 1:
        raise ValueError("At least one subset series is required")
    for executable in ("docker", "taskset", "lscpu", "rustc", "cargo", "/usr/bin/time"):
        if not shutil.which(executable):
            raise ValueError(f"Missing prerequisite: {executable}")
    time_version = capture(["/usr/bin/time", "--version"])
    if "GNU" not in time_version:
        raise ValueError("GNU /usr/bin/time is required")
    image = json.loads(capture(["docker", "image", "inspect", args.image]))[0]
    cpus = ",".join(map(str, available[:threads]))
    environment = {
        **THREAD_ENV, "RAYON_NUM_THREADS": str(threads), "OMP_NUM_THREADS": str(threads),
    }
    probe_tools = "panaroo uname python"
    blast_versions = ""
    if not args.no_concordance:
        probe_tools += " blastn makeblastdb"
        blast_versions = "blastn -version; makeblastdb -version; "
    probe = [
        "docker", "run", "--rm", "--pull=never", "--network=none",
        "--user", f"{os.getuid()}:{os.getgid()}", "--cpuset-cpus", cpus, image["Id"],
        "bash", "-euc",
        "if [ ! -x /usr/bin/time ]; then "
        "echo 'Missing prerequisite in the Panaroo image: /usr/bin/time (GNU time).' >&2; "
        "echo 'Build benchmarks/panaroo.Dockerfile and pass --image rustroo-panaroo-benchmark:1.5.1-time.' >&2; "
        "exit 127; fi; "
        f"for tool in {probe_tools}; do "
        "if ! command -v \"$tool\" >/dev/null; then "
        "echo \"Missing prerequisite in the Panaroo image: $tool\" >&2; exit 127; fi; done; "
        "/usr/bin/time --version; panaroo --version; "
        f"{blast_versions}uname -a; "
        "python -c 'import os; print(sorted(os.sched_getaffinity(0)))'",
    ]
    container_versions = capture(probe)
    if "GNU" not in container_versions:
        raise ValueError("The Panaroo image must provide GNU /usr/bin/time")
    root.mkdir(parents=True, exist_ok=False)
    directories = ["inputs", "manifests", "jobs", "protocol"]
    if not args.no_concordance:
        directories.append("concordance")
    for directory in directories:
        (root / directory).mkdir()
    shutil.copy2(binary, root / "rustroo")
    inputs = {}
    for name in names:
        destination = root / "inputs" / name
        shutil.copy2(source / name, destination)
        inputs[name] = {"sha256": digest(destination), "bytes": destination.stat().st_size}
    for key, selected in subsets.items():
        (root / "manifests" / f"{key}.txt").write_text("\n".join(selected) + "\n")
    script = Path(__file__).resolve()
    shutil.copy2(script, root / "protocol" / script.name)
    project = script.parent.parent
    source_hashes = {}
    for path in [project / "Cargo.toml", project / "Cargo.lock", *sorted((project / "src").rglob("*.rs"))]:
        destination = root / "protocol" / path.relative_to(project)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        source_hashes[str(path.relative_to(project))] = digest(destination)
    shutil.copy2(script.parent / "README.md", root / "protocol" / "README.md")
    shutil.copy2(script.parent / "panaroo.Dockerfile", root / "protocol" / "panaroo.Dockerfile")
    config = {
        "schema": 1, "run_id": hashlib.sha256(str(root).encode()).hexdigest()[:12],
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "sizes": sizes, "seeds": seeds, "repeats": args.repeats,
        "warmups": not args.no_warmups, "concordance": not args.no_concordance,
        "threads": threads, "cpus": cpus, "clean_mode": args.clean_mode,
        "image_tag": args.image, "image_id": image["Id"],
        "image_repo_digests": image.get("RepoDigests", []),
        "image_labels": image.get("Config", {}).get("Labels"),
        "binary_sha256": digest(root / "rustroo"), "driver_sha256": digest(script),
        "inputs": inputs, "subsets": subsets, "jobs": jobs,
        "tool_environment": environment,
        "blast": {"min_identity": 70.0, "min_coverage": 90.0, "evalue": "1e-10"},
    }
    save_json(root / "config.json", config)
    save_json(root / "environment.json", {
        "platform": platform.platform(), "python": sys.version,
        "rustroo_version": capture([str(root / "rustroo"), "--version"]),
        "rustc": capture(["rustc", "-Vv"]), "cargo": capture(["cargo", "-V"]),
        "docker": capture(["docker", "version"]), "cpu_topology": capture(["lscpu", "--json"]),
        "docker_daemon_resources": capture([
            "docker", "info", "--format",
            "NCPU={{.NCPU}}\nMemTotal={{.MemTotal}}\nOS={{.OperatingSystem}}\nKernel={{.KernelVersion}}",
        ]),
        "native_time": time_version, "container_versions": container_versions,
        "cpu_affinity_available": available, "source_sha256": source_hashes,
        "memory": Path("/proc/meminfo").read_text(),
        "build_environment": {key: os.environ.get(key) for key in ("RUSTFLAGS", "CARGO_PROFILE_RELEASE_LTO")},
    })
    measured = sum(not job["warmup"] for job in jobs)
    comparisons = 0 if args.no_concordance else len(subsets)
    print(f"Prepared {root}\n{measured} measured runs, {len(jobs) - measured} warmups, {comparisons} comparisons")


def verify_snapshot(root: Path, config: dict[str, Any], check_affinity: bool = True) -> None:
    if config["schema"] != 1 or digest(Path(__file__)) != config["driver_sha256"]:
        raise ValueError("Benchmark driver changed; use the saved protocol/benchmark.py")
    if digest(root / "rustroo") != config["binary_sha256"]:
        raise ValueError("Snapshotted binary changed")
    for name, information in config["inputs"].items():
        if digest(root / "inputs" / name) != information["sha256"]:
            raise ValueError(f"Snapshotted input changed: {name}")
    if check_affinity and not set(map(int, config["cpus"].split(","))).issubset(os.sched_getaffinity(0)):
        raise ValueError("The recorded CPU set is no longer available")


def run_job(root: Path, config: dict[str, Any], job: dict[str, Any]) -> None:
    folder = root / "jobs" / job["id"]
    if folder.exists():
        record = load_json(folder / "result.json")
        if record["status"] != "ok":
            raise ValueError(f"Failed attempt retained at {folder}; move it aside before retrying")
        current = output_metrics(folder / "output", job["tool"], config["subsets"][job["subset"]])
        if current != record["outputs"]:
            raise ValueError(f"Completed outputs changed: {folder}")
        print(f"Skipping completed {job['id']}", flush=True)
        return
    folder.mkdir()
    record: dict[str, Any] = {
        **job, "status": "running",
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host_load_before": Path("/proc/loadavg").read_text().strip(),
    }
    save_json(folder / "result.json", record)
    print(f"Running {job['id']}", flush=True)
    command = tool_command(config, root, job)
    container_name = f"rustroo_bench_{config['run_id']}_{job['id']}" if job["tool"] == "panaroo" else None
    try:
        record["launcher_wall_seconds"] = execute(
            command, folder, config["tool_environment"], container_name,
        )
        record["timing"] = timing_metrics(folder / "timing.txt")
        record["outputs"] = output_metrics(
            folder / "output", job["tool"], config["subsets"][job["subset"]],
        )
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError, KeyboardInterrupt) as error:
        record.update(status="failed", error=str(error) or type(error).__name__)
        save_json(folder / "result.json", record)
        raise
    record.update(
        status="ok", host_load_after=Path("/proc/loadavg").read_text().strip(),
        finished_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    )
    save_json(folder / "result.json", record)


def best_hits(path: Path) -> tuple[dict[str, list[str]], set[str], set[str]]:
    best: dict[str, tuple[float, dict[str, list[str]]]] = {}
    with path.open() as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != 13:
                raise ValueError(f"Expected 13 BLAST columns: {path}:{line_number}")
            qlen, slen = int(row[5]), int(row[6])
            qspan = abs(int(row[8]) - int(row[7])) + 1
            sspan = abs(int(row[10]) - int(row[9])) + 1
            if (
                qlen <= 0 or slen <= 0 or qspan > qlen or sspan > slen
                or not 1 <= int(row[7]) <= qlen or not 1 <= int(row[8]) <= qlen
                or not 1 <= int(row[9]) <= slen or not 1 <= int(row[10]) <= slen
                or not 0 <= int(row[4]) <= int(row[3]) or int(row[3]) <= 0
            ):
                raise ValueError(f"Invalid BLAST lengths: {path}:{line_number}")
            score = float(row[12])
            if not math.isfinite(score) or not 0 <= float(row[2]) <= 100:
                raise ValueError(f"Invalid BLAST score: {path}:{line_number}")
            query, subject = row[:2]
            previous = best.get(query)
            if previous is None or score > previous[0]:
                best[query] = (score, {subject: row})
            elif score == previous[0]:
                old = previous[1].get(subject)
                # Choose coverage deterministically when one subject has tied HSPs.
                def rank(hit: list[str]) -> tuple[float, float, tuple[str, ...]]:
                    coverage = min(
                        (abs(int(hit[8]) - int(hit[7])) + 1) / int(hit[5]),
                        (abs(int(hit[10]) - int(hit[9])) + 1) / int(hit[6]),
                    )
                    return coverage, float(hit[2]), tuple(hit)
                if old is None or rank(row) > rank(old):
                    previous[1][subject] = row
    ambiguous = {query for query, (_, hits) in best.items() if len(hits) > 1}
    unique = {
        query: next(iter(hits.values()))
        for query, (_, hits) in best.items() if len(hits) == 1
    }
    return unique, ambiguous, set(best)


def compare_references(
    forward_path: Path, reverse_path: Path,
    panaroo: dict[str, str], rustroo: dict[str, str],
    min_identity: float, min_coverage: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    forward, ambiguous_pan, hits_pan = best_hits(forward_path)
    reverse, ambiguous_rust, hits_rust = best_hits(reverse_path)
    if not hits_pan.issubset(panaroo) or not hits_rust.issubset(rustroo):
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
        qcoverage = 100 * (abs(int(row[8]) - int(row[7])) + 1) / int(row[5])
        scoverage = 100 * (abs(int(row[10]) - int(row[9])) + 1) / int(row[6])
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
        "accepted_pairs": len(pairs), "full_length_exact": sum(p["full_length_exact"] for p in pairs),
        "pairs_at_least_99pct_identity": sum(p["identity_pct"] >= 99 for p in pairs),
        "panaroo_accepted_fraction": len(pairs) / len(panaroo),
        "rustroo_accepted_fraction": len(pairs) / len(rustroo),
        "mean_identity_pct": statistics.mean(p["identity_pct"] for p in pairs) if pairs else None,
        "median_identity_pct": statistics.median(p["identity_pct"] for p in pairs) if pairs else None,
        "weighted_identity_pct": (
            100 * sum(p["identical_bases"] for p in pairs) / sum(p["alignment_columns"] for p in pairs)
        ) if pairs else None,
    }
    return metrics, pairs


def concordance(root: Path, config: dict[str, Any], subset: str) -> None:
    folder = root / "concordance" / subset
    for tool in TOOLS:
        job_folder = root / "jobs" / f"{subset}_r01_{tool}"
        result = load_json(job_folder / "result.json")
        if result["status"] != "ok":
            raise ValueError(f"Comparison requires successful repetition 1: {job_folder}")
        current = output_metrics(job_folder / "output", tool, config["subsets"][subset])
        if current != result["outputs"]:
            raise ValueError(f"Comparison source outputs changed: {job_folder}")
    pan = root / "jobs" / f"{subset}_r01_panaroo" / "output" / "pan_genome_reference.fa"
    rust = root / "jobs" / f"{subset}_r01_rustroo" / "output" / "genes.fa"
    source_hashes = {"panaroo": digest(pan), "rustroo": digest(rust)}
    if folder.exists():
        record = load_json(folder / "result.json")
        if record["status"] != "ok":
            raise ValueError(f"Failed comparison retained at {folder}; move it aside before retrying")
        if record["source_sha256"] != source_hashes:
            raise ValueError(f"Comparison reference inputs changed: {subset}")
        for name, expected in record["artifact_sha256"].items():
            if digest(folder / name) != expected:
                raise ValueError(f"Comparison artifact changed: {folder / name}")
        print(f"Skipping completed comparison {subset}", flush=True)
        return
    pan_sequences, rust_sequences = fasta_sequences(pan), fasta_sequences(rust)
    folder.mkdir()
    record: dict[str, Any] = {
        "status": "running", "subset": subset, "source_sha256": source_hashes,
    }
    save_json(folder / "result.json", record)
    script = """
set -euo pipefail
makeblastdb -in /work/panaroo.fa -dbtype nucl -out /work/panaroo_db
makeblastdb -in /work/rustroo.fa -dbtype nucl -out /work/rustroo_db
blastn -task blastn -query /work/panaroo.fa -db /work/rustroo_db -num_threads "$1" \
  -evalue "$2" -dust no -max_target_seqs "$3" -outfmt "$5" -out /work/forward.tsv
blastn -task blastn -query /work/rustroo.fa -db /work/panaroo_db -num_threads "$1" \
  -evalue "$2" -dust no -max_target_seqs "$4" -outfmt "$5" -out /work/reverse.tsv
"""
    shutil.copy2(pan, folder / "panaroo.fa")
    shutil.copy2(rust, folder / "rustroo.fa")
    container_name = f"rustroo_bench_{config['run_id']}_compare_{subset}"
    command = [
        *docker_base(config, root), "--name", container_name,
        "-v", f"{folder}:/work", "-w", "/work", config["image_id"],
        "bash", "-euc", script, "bash", str(config["threads"]), config["blast"]["evalue"],
        str(len(rust_sequences)), str(len(pan_sequences)), BLAST_FORMAT,
    ]
    try:
        execute(command, folder, config["tool_environment"], container_name)
        metrics, pairs = compare_references(
            folder / "forward.tsv", folder / "reverse.tsv", pan_sequences, rust_sequences,
            config["blast"]["min_identity"], config["blast"]["min_coverage"],
        )
        if pairs:
            write_tsv(folder / "shared_genes.tsv", pairs)
        else:
            (folder / "shared_genes.tsv").write_text(
                "panaroo_gene\trustroo_gene\tidentity_pct\tpanaroo_coverage_pct\t"
                "rustroo_coverage_pct\tidentical_bases\talignment_columns\tfull_length_exact\n"
            )
        record.update(
            status="ok", metrics=metrics,
            artifact_sha256={
                name: digest(folder / name)
                for name in ("forward.tsv", "reverse.tsv", "shared_genes.tsv")
            },
        )
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError, KeyboardInterrupt) as error:
        record.update(status="failed", error=str(error) or type(error).__name__)
        save_json(folder / "result.json", record)
        raise
    save_json(folder / "result.json", record)


def spread(values: list[float]) -> dict[str, float]:
    return {"median": statistics.median(values), "min": min(values), "max": max(values)}


def svg_plot(
    path: Path, title: str, ylabel: str, series: dict[str, list[dict[str, float]]],
) -> None:
    points = [point for values in series.values() for point in values]
    sizes = sorted({point["n"] for point in points})
    maximum = max(point["max"] for point in points) * 1.1 or 1
    left, top, width, height = 90, 55, 660, 300
    x = lambda n: left + (n - min(sizes)) / (max(sizes) - min(sizes) or 1) * width
    y = lambda value: top + height * (1 - value / maximum)
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="430" viewBox="0 0 900 430">',
        '<rect width="900" height="430" fill="white"/>',
        '<g font-family="sans-serif" font-size="13">',
        f'<text x="90" y="28" font-size="18">{html.escape(title)}</text>',
        f'<text x="90" y="420">{html.escape(ylabel)}; bars: range across subset medians</text>',
    ]
    for tick in range(6):
        value = maximum * tick / 5
        lines += [
            f'<line x1="{left}" y1="{y(value)}" x2="{left + width}" y2="{y(value)}" stroke="#ddd"/>',
            f'<text x="82" y="{y(value) + 4}" text-anchor="end">{value:.3g}</text>',
        ]
    for n in sizes:
        lines.append(f'<text x="{x(n)}" y="377" text-anchor="middle">{n}</text>')
    lines.append('<text x="420" y="398" text-anchor="middle">Chromosome GFF files</text>')
    for index, (name, values) in enumerate(series.items()):
        color = ("#0072B2", "#D55E00", "#009E73")[index % 3]
        coords = " ".join(f"{x(point['n'])},{y(point['median'])}" for point in values)
        lines.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2"/>')
        for point in values:
            px = x(point["n"])
            lines.append(f'<line x1="{px}" x2="{px}" y1="{y(point["min"])}" y2="{y(point["max"])}" stroke="{color}"/>')
            lines.append(f'<circle cx="{px}" cy="{y(point["median"])}" r="4" fill="{color}"/>')
        lines.append(f'<text x="770" y="{80 + index * 22}" fill="{color}">{html.escape(name)}</text>')
    lines.append("</g></svg>\n")
    path.write_text("\n".join(lines))


def report(root: Path, config: dict[str, Any]) -> None:
    include_concordance = config.get("concordance", True)
    raw = []
    for job in config["jobs"]:
        record = load_json(root / "jobs" / job["id"] / "result.json")
        if record["status"] != "ok":
            raise ValueError(f"Incomplete benchmark: {job['id']}")
        if not job["warmup"]:
            raw.append({
                **job, **record["timing"], **record["outputs"],
                "launcher_wall_seconds": record["launcher_wall_seconds"],
                "started_utc": record["started_utc"],
            })
    output = root / "report"
    output.mkdir(exist_ok=True)
    write_tsv(output / "runs.tsv", raw)
    subset_rows = []
    for subset in sorted(config["subsets"]):
        for tool in TOOLS:
            rows = [row for row in raw if row["subset"] == subset and row["tool"] == tool]
            if len(rows) != config["repeats"]:
                raise ValueError(f"Wrong repeat count: {subset} {tool}")
            summary = {
                "subset": subset, "n": rows[0]["n"], "series": rows[0]["series"],
                "seed": rows[0]["seed"], "tool": tool, "repeats": len(rows),
                "reference_hash_variants": len({row["reference_sha256"] for row in rows}),
                "presence_hash_variants": len({row["presence_file_sha256"] for row in rows}),
            }
            for metric in ("wall_seconds", "max_process_rss_kib", "average_cpu_cores", "groups", "references"):
                summary.update({f"{metric}_{key}": value for key, value in spread([row[metric] for row in rows]).items()})
            subset_rows.append(summary)
    write_tsv(output / "subset_summary.tsv", subset_rows)
    aggregate = []
    speedups = []
    concordances = []
    for n in config["sizes"]:
        for tool in TOOLS:
            selected = [row for row in subset_rows if row["n"] == n and row["tool"] == tool]
            row = {"n": n, "tool": tool, "subset_series": len(selected)}
            for metric in ("wall_seconds", "max_process_rss_kib", "average_cpu_cores", "groups", "references"):
                row.update({
                    f"{metric}_{key}": value
                    for key, value in spread([r[f"{metric}_median"] for r in selected]).items()
                })
            aggregate.append(row)
        ratios = []
        for subset in sorted(config["subsets"]):
            selected = {r["tool"]: r for r in subset_rows if r["subset"] == subset and r["n"] == n}
            if selected:
                ratio = selected["panaroo"]["wall_seconds_median"] / selected["rustroo"]["wall_seconds_median"]
                ratios.append(ratio)
                speedups.append({"subset": subset, "n": n, "panaroo_over_rustroo": ratio})
                if include_concordance:
                    record = load_json(root / "concordance" / subset / "result.json")
                    if record["status"] != "ok":
                        raise ValueError(f"Incomplete concordance: {subset}")
                    concordances.append({"subset": subset, "n": n, **record["metrics"]})
        print(f"{n:2d} genomes: median paired speedup {statistics.median(ratios):.3f}x")
    write_tsv(output / "scaling_summary.tsv", aggregate)
    write_tsv(output / "paired_speedups.tsv", speedups)
    write_tsv(output / "speedup_summary.tsv", [
        {
            "n": n, "subset_series": len(config["seeds"]),
            **spread([row["panaroo_over_rustroo"] for row in speedups if row["n"] == n]),
        } for n in config["sizes"]
    ])
    if include_concordance:
        write_tsv(output / "concordance.tsv", concordances)
        concordance_summary = []
        for n in config["sizes"]:
            selected = [row for row in concordances if row["n"] == n]
            row = {"n": n, "subset_series": len(selected)}
            for metric in (
                "accepted_pairs", "full_length_exact", "panaroo_accepted_fraction",
                "rustroo_accepted_fraction", "mean_identity_pct", "weighted_identity_pct",
            ):
                values = [r[metric] for r in selected if r[metric] is not None]
                row[f"{metric}_subsets_with_value"] = len(values)
                for key in ("median", "min", "max"):
                    row[f"{metric}_{key}"] = spread(values)[key] if values else None
            concordance_summary.append(row)
        write_tsv(output / "concordance_summary.tsv", concordance_summary)
    for filename, title, metric, units in (
        ("runtime.svg", "Runtime scaling", "wall_seconds", "Elapsed seconds"),
        ("max_process_rss.svg", "Maximum-process RSS (not aggregate pipeline RAM)", "max_process_rss_kib", "KiB"),
        ("groups.svg", "Gene-family counts", "groups", "Groups"),
    ):
        svg_plot(output / filename, title, units, {
            tool: [
                {"n": row["n"], **{key: row[f"{metric}_{key}"] for key in ("median", "min", "max")}}
                for row in aggregate if row["tool"] == tool
            ] for tool in TOOLS
        })
    svg_plot(output / "speedup.svg", "Paired runtime ratio", "Panaroo / rustroo", {
        "speedup": [
            {"n": n, **spread([row["panaroo_over_rustroo"] for row in speedups if row["n"] == n])}
            for n in config["sizes"]
        ],
    })
    if include_concordance:
        svg_plot(output / "matched_references.svg", "Accepted reciprocal reference matches", "Fraction of references", {
            tool: [
                {"n": n, **spread([row[f"{tool}_accepted_fraction"] for row in concordances if row["n"] == n])}
                for n in config["sizes"]
            ] for tool in TOOLS
        })
    changes = [r for r in subset_rows if r["reference_hash_variants"] > 1 or r["presence_hash_variants"] > 1]
    save_json(output / "output_variability.json", {
        "variable_subsets": changes,
        "concordance_uses_repeat": 1 if include_concordance else None,
        "note": "Reference hashes normalize FASTA order; presence hashes are byte-level.",
    })
    if changes:
        print(f"WARNING: output hashes vary in {len(changes)} subset/tool combinations; inspect output_variability.json")
    print(f"Tables and SVG figures: {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    create = commands.add_parser("prepare", help="Snapshot inputs, binary, versions and randomized schedule")
    create.add_argument("--input", type=Path, default=Path("gff_test/chromosome"))
    create.add_argument("--binary", type=Path, default=Path("target/release/rustroo"))
    create.add_argument("--output", type=Path, required=True)
    create.add_argument("--image", default="staphb/panaroo:1.5.1")
    create.add_argument("--threads", type=int)
    create.add_argument("--seed", type=int, default=1729)
    create.add_argument("--subsets", type=int, default=5)
    create.add_argument("--repeats", type=int, default=5)
    create.add_argument("--sizes", type=int, nargs="+", default=list(SIZES))
    create.add_argument("--no-warmups", action="store_true")
    create.add_argument("--no-concordance", action="store_true")
    create.add_argument("--clean-mode", choices=("strict", "moderate", "sensitive"), default="strict")
    for name in ("run", "report"):
        command = commands.add_parser(name)
        command.add_argument("--output", type=Path, required=True)
        if name == "run":
            command.add_argument("--stage", choices=("performance", "concordance", "all"), default="all")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args)
        return
    root = args.output.resolve()
    config = load_json(root / "config.json")
    if args.action == "run":
        if sys.platform != "linux":
            raise ValueError("Execute the benchmark in Linux/WSL")
        if args.stage == "concordance" and not config.get("concordance", True):
            raise ValueError("Reference concordance was disabled when this experiment was prepared")
        import fcntl

        with (root / "execution.lock").open("a") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise ValueError("Another benchmark stage is running in this output directory") from error
            verify_snapshot(root, config)
            if args.stage in ("performance", "all"):
                for job in config["jobs"]:
                    run_job(root, config, job)
            if args.stage in ("concordance", "all") and config.get("concordance", True):
                for subset in sorted(config["subsets"]):
                    concordance(root, config, subset)
            if args.stage == "all":
                report(root, config)
        return
    verify_snapshot(root, config, check_affinity=False)
    report(root, config)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"benchmark: {error}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("benchmark: interrupted; failed attempt retained, no results imputed", file=sys.stderr)
        sys.exit(130)
