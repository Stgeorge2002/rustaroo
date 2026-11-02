#!/usr/bin/env python3
"""UniProt length accuracy across the full Anatra ST6 scaling grid.

Design:

1. For each species, read the n=3000 rustroo and Panaroo reference FASTAs.
2. Take gene names present in both (shared n=3000 set) — one unique list per
   species.
3. Look those names up once in UniProt (per-species taxonomy id; disk-cached).
4. Score every available cohort size (10 … 3000) for both tools against that
   same UniProt length table: a gene in a smaller run is scored only when it
   is in the n=3000 shared set and UniProt knows the name.
5. Emit one cross-species table so rustroo vs Panaroo can be compared across
   all sizes at once.

Match rule: reference bp equals UniProt amino acids × 3, or that plus the
stop codon. Names UniProt does not return are excluded from the match rate.

Example:

  python3 benchmarks/uniprot_lengths_st6.py \\
    --results-root anatra_st6_results/st6 \\
    --output anatra_st6_results/uniprot_lengths
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
TRANSIENT_HTTP = {429, 500, 502, 503, 504}

# NCBI taxonomy ids for gene:"name" AND taxonomy_id:N UniProt queries.
SPECIES_TAXON = {
    "campylobacter_jejuni": 197,
    "streptococcus_pneumoniae": 1313,
    "staphylococcus_aureus": 1280,
    "escherichia_coli": 562,
    "klebsiella_pneumoniae": 573,
    "pseudomonas_aeruginosa": 287,
}

# Matches anatra_st6_bench.sh SIZE case.
SIZES = (10, 25, 50, 100, 200, 350, 500, 750, 1015, 2000, 3000)
LOOKUP_SIZE = 3000

TOOLS = ("rustroo", "panaroo")
REF_NAME = {
    "rustroo": "genes.fa",
    "panaroo": "pan_genome_reference.fa",
}


def read_lengths(path: Path):
    name = None
    length = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        yield name, length


def request(url: str, retries: int = 12) -> str:
    last_error: BaseException | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "rustroo-uniprot-lengths-st6/1.0"}
            )
            with urllib.request.urlopen(req, timeout=60) as response:
                return response.read().decode()
        except urllib.error.HTTPError as error:
            body = error.read().decode(errors="replace")[:300]
            last_error = error
            if error.code not in TRANSIENT_HTTP or attempt + 1 == retries:
                raise RuntimeError("HTTP %s: %s" % (error.code, body)) from error
            retry_after = error.headers.get("Retry-After")
            try:
                wait = float(retry_after) if retry_after else min(2 ** attempt, 60)
            except ValueError:
                wait = min(2 ** attempt, 60)
            print(
                "  retry %s/%s HTTP %s; sleep %.0fs"
                % (attempt + 1, retries, error.code, wait),
                flush=True,
            )
            time.sleep(wait)
            continue
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            last_error = error
            if attempt + 1 == retries:
                raise
            wait = min(2 ** attempt, 60)
            print(
                "  retry %s/%s %s; sleep %.0fs"
                % (attempt + 1, retries, type(error).__name__, wait),
                flush=True,
            )
            time.sleep(wait)
            continue
    raise RuntimeError(url) from last_error


def lookup(name: str, taxon: int, retries: int = 12) -> dict:
    quoted = name.replace('"', "")
    query = '(gene:"%s") AND taxonomy_id:%s' % (quoted, taxon)
    encoded = urllib.parse.urlencode({
        "query": query,
        "fields": "accession,length,gene_names",
        "format": "tsv",
        "size": "25",
    })
    lines = request(SEARCH + "?" + encoded, retries=retries).splitlines()[1:]
    lengths = []
    accessions = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2 or not parts[1].isdigit():
            continue
        accessions.append(parts[0])
        lengths.append(int(parts[1]))
    if not lengths:
        return {"status": "absent"}
    modal = Counter(lengths).most_common(1)[0][0]
    return {
        "status": "found",
        "uniprot": accessions[0],
        "uniprot_aa": modal,
        "n_entries": len(lengths),
        "n_lengths": len(set(lengths)),
    }


def agrees(reference_bp: int, amino_acids: int) -> bool:
    return reference_bp in (amino_acids * 3, amino_acids * 3 + 3)


def safe_cache_name(name: str) -> str:
    slug = "".join(character if character.isalnum() else "_" for character in name)[:80]
    digest = hashlib.sha256(name.encode("utf-8")).hexdigest()[:12]
    return "%s_%s" % (slug, digest)


def find_ref(run_dir: Path, tool: str) -> Path:
    path = run_dir / REF_NAME[tool]
    if path.is_file():
        return path
    nested = run_dir / "out" / REF_NAME[tool]
    if nested.is_file():
        return nested
    raise FileNotFoundError(
        "missing %s under %s (checked run dir and out/)" % (REF_NAME[tool], run_dir)
    )


def run_dir(results_root: Path, species: str, size: int, tool: str, rep: int) -> Path:
    return (
        results_root
        / species
        / ("n%04d" % size)
        / ("%s_rep%02d" % (tool, rep))
    )


def discover_species(results_root: Path, wanted: list[str] | None) -> list[str]:
    found = []
    for species in SPECIES_TAXON:
        if wanted is not None and species not in wanted:
            continue
        if (results_root / species).is_dir():
            found.append(species)
    if wanted is not None:
        missing = [name for name in wanted if name not in found]
        if missing:
            raise SystemExit("species dirs not under %s: %s" % (results_root, ", ".join(missing)))
    if not found:
        raise SystemExit("no species directories under %s" % results_root)
    return found


def load_tool_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for name, bp in read_lengths(path):
        lengths.setdefault(name, bp)
    return lengths


def rate(bucket: Counter) -> float:
    known = bucket["found"]
    return 100.0 * bucket["match"] / known if known else 0.0


def shared_n3000_names(
    results_root: Path, species: str, rep: int
) -> tuple[list[str], dict[str, dict[str, int]]]:
    by_tool: dict[str, dict[str, int]] = {}
    for tool in TOOLS:
        path = find_ref(run_dir(results_root, species, LOOKUP_SIZE, tool, rep), tool)
        by_tool[tool] = load_tool_lengths(path)
    names = sorted(set(by_tool["rustroo"]) & set(by_tool["panaroo"]))
    return names, by_tool


def fetch_uniprot(
    species: str,
    taxon: int,
    names: list[str],
    cache_dir: Path,
    retries: int,
) -> tuple[dict[str, dict], int]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    found: dict[str, dict] = {}
    pending = 0
    for index, name in enumerate(names, 1):
        path = cache_dir / (safe_cache_name(name) + ".json")
        if path.is_file():
            row = json.loads(path.read_text())
        else:
            while True:
                try:
                    row = lookup(name, taxon, retries=retries)
                    break
                except RuntimeError as error:
                    print(
                        "  UniProt still failing on %s (%s); wait 90s then retry"
                        % (name, error),
                        flush=True,
                    )
                    time.sleep(90)
            path.write_text(json.dumps(row))
            pending += 1
            time.sleep(0.34)
        found[name] = row
        print(
            "  lookup %s/%s %s %s"
            % (index, len(names), name, row.get("status", "absent")),
            flush=True,
        )
    return found, pending


def score_run(
    species: str,
    size: int,
    tool: str,
    lengths: dict[str, int],
    shared_names: set[str],
    uniprot: dict[str, dict],
) -> tuple[list[list[str]], Counter, Counter]:
    """Return detail rows, per-tool counts, and shared-name contingency counts."""
    rows: list[list[str]] = []
    counts = Counter()
    contingency = Counter()
    for name in sorted(set(lengths) & shared_names):
        row = uniprot[name]
        amino = row.get("uniprot_aa")
        status = row.get("status", "absent")
        bp = lengths[name]
        match = status == "found" and agrees(bp, amino)
        counts["references"] += 1
        counts[status] += 1
        counts["match"] += int(match)
        rows.append([
            species,
            str(size),
            tool,
            name,
            str(bp),
            status,
            row.get("uniprot", ""),
            str(amino or ""),
            str(amino * 3 + 3 if amino else ""),
            str(row.get("n_entries", "")),
            str(row.get("n_lengths", "")),
            "yes" if match else "no",
        ])
    # Contingency needs both tools; filled by caller for a size.
    return rows, counts, contingency


def contingency_for_size(
    rust_lengths: dict[str, int],
    pan_lengths: dict[str, int],
    shared_names: set[str],
    uniprot: dict[str, dict],
) -> Counter:
    out = Counter()
    for name in shared_names:
        if name not in rust_lengths or name not in pan_lengths:
            continue
        row = uniprot[name]
        if row.get("status") != "found":
            continue
        amino = row["uniprot_aa"]
        r_ok = agrees(rust_lengths[name], amino)
        p_ok = agrees(pan_lengths[name], amino)
        if r_ok and p_ok:
            out["both_match"] += 1
        elif r_ok:
            out["rustroo_only"] += 1
        elif p_ok:
            out["panaroo_only"] += 1
        else:
            out["neither"] += 1
    return out


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root",
        type=Path,
        default=repo / "anatra_st6_results" / "st6",
        help="ST6 tree with <species>/nNNNN/<tool>_repNN/",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repo / "anatra_st6_results" / "uniprot_lengths",
    )
    parser.add_argument(
        "--sizes",
        type=int,
        nargs="+",
        default=list(SIZES),
        help="cohort sizes to score (default: full ST6 grid)",
    )
    parser.add_argument(
        "--lookup-size",
        type=int,
        default=LOOKUP_SIZE,
        help="cohort that defines the shared gene list (default: 3000)",
    )
    parser.add_argument("--rep", type=int, default=1, help="replicate number")
    parser.add_argument(
        "--species",
        nargs="+",
        choices=sorted(SPECIES_TAXON),
        default=None,
        help="subset of species; default = every species dir present",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="cap shared n=3000 names looked up per species; 0 means all",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=12,
        help="transient HTTP retries per UniProt request before a longer pause",
    )
    args = parser.parse_args()

    if not args.results_root.is_dir():
        raise SystemExit("results root not found: %s" % args.results_root)
    if args.lookup_size not in args.sizes:
        # Still allow scoring a subset, but the gene list always comes from
        # --lookup-size, which must exist on disk.
        pass

    species_list = discover_species(args.results_root, args.species)
    args.output.mkdir(parents=True, exist_ok=True)

    detail_header = [
        "species", "size", "tool", "gene", "reference_bp", "status", "uniprot",
        "uniprot_aa", "uniprot_bp", "n_entries", "n_lengths", "length_match",
    ]
    summary_header = [
        "species", "size", "tool", "genes_scored", "uniprot_found",
        "length_matches", "match_pct", "both_match", "rustroo_only",
        "panaroo_only", "neither",
    ]
    gene_header = [
        "species", "gene", "status", "uniprot", "uniprot_aa", "uniprot_bp",
        "n_entries", "n_lengths", "rustroo_n3000_bp", "panaroo_n3000_bp",
    ]

    all_detail: list[list[str]] = []
    all_summary: list[list[str]] = []
    all_genes: list[list[str]] = []
    overall = {tool: Counter() for tool in TOOLS}
    overall_by_size = {(size, tool): Counter() for size in args.sizes for tool in TOOLS}

    for species in species_list:
        taxon = SPECIES_TAXON[species]
        print("=" * 72, flush=True)
        print("species=%s taxon=%s lookup_size=%s" % (species, taxon, args.lookup_size), flush=True)

        names, n3000 = shared_n3000_names(args.results_root, species, args.rep)
        if args.limit:
            names = names[: args.limit]
        shared_set = set(names)
        print(
            "  shared n=%s genes=%s (rustroo=%s panaroo=%s)"
            % (
                args.lookup_size,
                len(names),
                len(n3000["rustroo"]),
                len(n3000["panaroo"]),
            ),
            flush=True,
        )

        uniprot, pending = fetch_uniprot(
            species=species,
            taxon=taxon,
            names=names,
            cache_dir=args.output / species / "cache",
            retries=args.retries,
        )
        print("  new UniProt lookups=%s" % pending, flush=True)

        for name in names:
            row = uniprot[name]
            amino = row.get("uniprot_aa")
            all_genes.append([
                species,
                name,
                row.get("status", ""),
                row.get("uniprot", ""),
                str(amino or ""),
                str(amino * 3 + 3 if amino else ""),
                str(row.get("n_entries", "")),
                str(row.get("n_lengths", "")),
                str(n3000["rustroo"].get(name, "")),
                str(n3000["panaroo"].get(name, "")),
            ])

        species_detail: list[list[str]] = []
        species_summary: list[list[str]] = []

        for size in args.sizes:
            by_tool: dict[str, dict[str, int]] = {}
            missing = False
            for tool in TOOLS:
                directory = run_dir(args.results_root, species, size, tool, args.rep)
                if not directory.is_dir():
                    print("  skip missing %s" % directory, flush=True)
                    missing = True
                    break
                try:
                    by_tool[tool] = load_tool_lengths(find_ref(directory, tool))
                except FileNotFoundError as error:
                    print("  skip %s" % error, flush=True)
                    missing = True
                    break
            if missing:
                continue

            size_counts = {}
            for tool in TOOLS:
                rows, counts, _ = score_run(
                    species=species,
                    size=size,
                    tool=tool,
                    lengths=by_tool[tool],
                    shared_names=shared_set,
                    uniprot=uniprot,
                )
                species_detail.extend(rows)
                all_detail.extend(rows)
                size_counts[tool] = counts
                overall[tool].update(counts)
                overall_by_size[(size, tool)].update(counts)

            contingency = contingency_for_size(
                by_tool["rustroo"], by_tool["panaroo"], shared_set, uniprot
            )
            for tool in TOOLS:
                bucket = size_counts[tool]
                species_summary.append([
                    species,
                    str(size),
                    tool,
                    str(bucket["references"]),
                    str(bucket["found"]),
                    str(bucket["match"]),
                    "%.4f" % rate(bucket),
                    str(contingency["both_match"]),
                    str(contingency["rustroo_only"]),
                    str(contingency["panaroo_only"]),
                    str(contingency["neither"]),
                ])
                all_summary.append(species_summary[-1])

            print(
                "  n=%s  rustroo %.1f%% (%s/%s)  panaroo %.1f%% (%s/%s)  both/r/p/neither=%s/%s/%s/%s"
                % (
                    size,
                    rate(size_counts["rustroo"]),
                    size_counts["rustroo"]["match"],
                    size_counts["rustroo"]["found"],
                    rate(size_counts["panaroo"]),
                    size_counts["panaroo"]["match"],
                    size_counts["panaroo"]["found"],
                    contingency["both_match"],
                    contingency["rustroo_only"],
                    contingency["panaroo_only"],
                    contingency["neither"],
                ),
                flush=True,
            )

        write_tsv(args.output / species / "lengths.tsv", detail_header, species_detail)
        write_tsv(args.output / species / "summary.tsv", summary_header, species_summary)
        write_tsv(
            args.output / species / "shared_genes.tsv",
            gene_header,
            [row for row in all_genes if row[0] == species],
        )

    write_tsv(args.output / "lengths.tsv", detail_header, all_detail)
    write_tsv(args.output / "summary.tsv", summary_header, all_summary)
    write_tsv(args.output / "shared_genes.tsv", gene_header, all_genes)

    by_size_rows = []
    for size in args.sizes:
        for tool in TOOLS:
            bucket = overall_by_size[(size, tool)]
            by_size_rows.append([
                str(size),
                tool,
                str(bucket["references"]),
                str(bucket["found"]),
                str(bucket["match"]),
                "%.4f" % rate(bucket),
            ])
    write_tsv(
        args.output / "by_size.tsv",
        ["size", "tool", "genes_scored", "uniprot_found", "length_matches", "match_pct"],
        by_size_rows,
    )

    print()
    print("wrote %s" % (args.output / "summary.tsv"))
    print("wrote %s" % (args.output / "by_size.tsv"))
    print("overall across all species × sizes:")
    for tool in TOOLS:
        bucket = overall[tool]
        print(
            "  %s  %s/%s (%.1f%%) of UniProt-known shared genes"
            % (tool, bucket["match"], bucket["found"], rate(bucket))
        )
    print("per size:")
    for size in args.sizes:
        rust = overall_by_size[(size, "rustroo")]
        pan = overall_by_size[(size, "panaroo")]
        if not rust["references"] and not pan["references"]:
            continue
        print(
            "  n=%-4s  rustroo %.1f%% (%s/%s)  panaroo %.1f%% (%s/%s)"
            % (
                size,
                rate(rust),
                rust["match"],
                rust["found"],
                rate(pan),
                pan["match"],
                pan["found"],
            )
        )


if __name__ == "__main__":
    main()
