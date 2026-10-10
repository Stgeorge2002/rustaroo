#!/usr/bin/env python3
"""Annotation-level family concordance, separate from timed pangenome builds."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from typing import Iterable, NamedTuple, Optional

from benchmark import BLAST_FORMAT, digest, fasta_sequences, save_json


Key = tuple[str, str]
FIELDS = {
    "components.tsv": (
        "component", "relationship", "rustroo_groups", "panaroo_groups", "group_delta",
        "shared_annotations", "rustroo_cooccurring_samples", "rustroo_multilocus_samples",
        "sequence_supported_membership_pairs", "membership_pairs", "review_evidence",
    ),
    "families.tsv": (
        "tool", "gene", "component", "relationship", "samples", "source_annotations",
        "excluded_source_annotations", "shared_annotations", "shared_fraction", "unshared_annotations", "refound_members",
        "fusion_pieces", "missing_id_members", "source_sequence_mismatches",
        "source_contig_mismatches", "reference_bp", "top_sequence_hits", "top_hit_tied",
    ),
    "membership_pairs.tsv": (
        "rustroo_gene", "panaroo_gene", "shared_annotations", "rustroo_shared_fraction",
        "panaroo_shared_fraction", "membership_jaccard", "presence_jaccard",
        "source_sequence_exact", "source_sequence_mismatches", "source_contig_mismatches",
        "reference_status", "identity_pct", "rustroo_coverage_pct", "panaroo_coverage_pct",
        "passes_sequence_thresholds",
    ),
    "sequence_hits.tsv": (
        "rustroo_gene", "panaroo_gene", "identity_pct", "rustroo_coverage_pct",
        "panaroo_coverage_pct", "bitscore", "passes_sequence_thresholds",
        "is_top_bitscore", "has_shared_annotations",
    ),
    "annotation_exceptions.tsv": ("tool", "gene", "sample", "annotation_id", "reason"),
}


@dataclass(frozen=True)
class Member:
    gene: str
    sample: str
    annotation: str
    contig: str
    start: int
    end: int
    strand: int
    locus: str
    origin: str
    allele: str


class PanSource(NamedTuple):
    contig: str
    dna: str


@dataclass(frozen=True)
class Hit:
    identity: float
    rust_coverage: float
    pan_coverage: float
    bitscore: float

    def passes(self, identity: float, coverage: float) -> bool:
        return self.identity >= identity and min(self.rust_coverage, self.pan_coverage) >= coverage


def progress(message: str) -> None:
    print(f"family audit: {message}", file=sys.stderr, flush=True)


def table(path: Path, fields: Iterable[str], delimiter: str = "\t") -> Iterable[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        header = reader.fieldnames
        if header is None or len(header) != len(set(header)) or not set(fields).issubset(header):
            raise ValueError(f"Missing or duplicate columns in {path}; need {sorted(fields)}")
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed row: {path}:{number}")
            yield row


def load_presence(path: Path) -> tuple[list[str], dict[str, set[str]]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        samples = header[1:]
        if (
            len(header) < 2 or header[0] != "gene" or not all(samples)
            or len(samples) != len(set(samples))
        ):
            raise ValueError(f"Invalid presence header: {path}")
        presence = {}
        for number, row in enumerate(reader, 2):
            if len(row) != len(header) or not row[0] or row[0] in presence:
                raise ValueError(f"Invalid/duplicate presence row: {path}:{number}")
            if any(value not in {"0", "1"} for value in row[1:]):
                raise ValueError(f"Nonbinary presence row: {path}:{number}")
            presence[row[0]] = {sample for sample, value in zip(samples, row[1:]) if value == "1"}
        if not presence or any(not values for values in presence.values()):
            raise ValueError(f"Empty family or matrix: {path}")
    return samples, presence


def load_members(
    path: Path, presence: dict[str, set[str]], alleles: dict[str, str],
) -> list[Member]:
    members = []
    observed: dict[str, set[str]] = defaultdict(set)
    loci: dict[str, tuple[str, str, str, int]] = {}
    seen = set()
    for row in table(path, ("gene", "sample", "annotation_id", "contig", "start", "end",
                            "strand", "locus", "origin", "allele")):
        gene, sample = row["gene"], row["sample"]
        if gene not in presence or sample not in presence[gene]:
            raise ValueError(f"Membership does not match presence.tsv: {gene}/{sample}")
        if row["allele"] not in alleles:
            raise ValueError(f"Unknown member allele: {row['allele']}")
        member = Member(
            gene, sample, row["annotation_id"], row["contig"], int(row["start"]), int(row["end"]),
            int(row["strand"]), row["locus"], row["origin"], row["allele"],
        )
        if (
            member.start < 1 or member.end < member.start or member.strand not in {-1, 1}
            or not member.contig or not member.locus
            or member.origin not in {"annotated", "fusion_piece", "refound"}
            or (member.origin == "refound" and member.annotation)
        ):
            raise ValueError(f"Invalid member: {member}")
        if member in seen:
            raise ValueError(f"Duplicate member: {member}")
        seen.add(member)
        location = (sample, member.contig, gene, member.strand)
        if member.locus in loci and loci[member.locus] != location:
            raise ValueError(f"Physical locus reused across incompatible locations: {member.locus}")
        loci[member.locus] = location
        observed[gene].add(sample)
        members.append(member)
    if dict(observed) != presence:
        raise ValueError("Membership sample sets differ from presence.tsv")
    return members


def load_pan(
    folder: Path, samples: list[str],
) -> tuple[dict[str, set[Key]], dict[str, set[str]], dict[Key, PanSource], list[dict[str, str]]]:
    matrix = folder / "gene_presence_absence.csv"
    with matrix.open(newline="") as handle:
        header = next(csv.reader(handle), [])
    metadata = {"Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
                "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
                "Accessory Fragment", "Accessory Order with Fragment", "QC",
                "Min group size nuc", "Max group size nuc", "Avg group size nuc"}
    if set(header) - metadata != set(samples):
        raise ValueError("Panaroo and Rustroo matrix sample sets differ; use identical GFF cohorts")
    progress("Reading Panaroo gene_data.csv (pooled exact DNA sequences)")
    data: dict[Key, PanSource] = {}
    sample_names = {sample: sample for sample in samples}
    dna_pool: dict[str, str] = {}
    contig_pool: dict[str, str] = {}
    for count, row in enumerate(table(folder / "gene_data.csv",
                                     ("gff_file", "annotation_id", "scaffold_name", "dna_sequence"), ","), 1):
        sample = sample_names.get(row["gff_file"])
        key = (row["gff_file"], row["annotation_id"])
        if sample is None or not key[1] or key in data:
            raise ValueError(f"Unexpected sample, missing ID or duplicate gene_data annotation: {key}")
        if not row["dna_sequence"] or not row["scaffold_name"]:
            raise ValueError(f"Missing sequence/scaffold for {key}")
        dna = row["dna_sequence"].upper()
        contig = row["scaffold_name"]
        # Share identical strings without substituting hash equality for DNA equality.
        data[sample, key[1]] = PanSource(
            contig_pool.setdefault(contig, contig), dna_pool.setdefault(dna, dna),
        )
        if count % 50_000 == 0:
            progress(f"Panaroo gene data: {count:,} annotations, {len(dna_pool):,} distinct DNA sequences")
    progress(f"Panaroo gene data complete: {len(data):,} annotations, {len(dna_pool):,} distinct DNA sequences")
    progress("Reading Panaroo family memberships")
    groups: dict[str, set[Key]] = {}
    presence: dict[str, set[str]] = {}
    exceptions = []
    for row in table(matrix, ("Gene", *samples), ","):
        gene = row["Gene"]
        if not gene or gene in groups:
            raise ValueError(f"Missing/duplicate Panaroo group: {gene}")
        groups[gene] = set()
        presence[gene] = set()
        for sample in samples:
            tokens = row[sample].split(";") if row[sample] else []
            for token in tokens:
                token = token.strip()
                if not token:
                    raise ValueError(f"Empty annotation token in {gene}/{sample}")
                key = (sample, token)
                presence[gene].add(sample)
                if "_refound_" in token:
                    exceptions.append({"tool": "panaroo", "gene": gene, "sample": sample,
                                       "annotation_id": token, "reason": "refound"})
                elif key not in data:
                    raise ValueError(f"Panaroo matrix ID absent from gene_data.csv: {key}")
                elif key in groups[gene]:
                    raise ValueError(f"Repeated Panaroo annotation within group {gene}: {key}")
                else:
                    groups[gene].add(key)
        if not presence[gene]:
            raise ValueError(f"Empty Panaroo group: {gene}")
    return groups, presence, data, exceptions


def load_hits(path: Path, rust: dict[str, str], pan: dict[str, str]) -> dict[tuple[str, str], Hit]:
    hits = {}
    with path.open() as handle:
        for number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != 13:
                raise ValueError(f"Expected {BLAST_FORMAT}: {path}:{number}")
            r, p = row[:2]
            if r not in rust or p not in pan:
                raise ValueError(f"BLAST identifiers do not match references: {path}:{number}")
            identity, bitscore, evalue = float(row[2]), float(row[12]), float(row[11])
            length, nident, qlen, slen = map(int, row[3:7])
            qs, qe, ss, se = map(int, row[7:11])
            if (
                qlen != len(rust[r]) or slen != len(pan[p]) or not 0 <= identity <= 100
                or not all(math.isfinite(value) for value in (identity, bitscore, evalue))
                or bitscore < 0 or evalue < 0 or not 0 <= nident <= length or length <= 0
                or not 1 <= qs <= qlen or not 1 <= qe <= qlen
                or not 1 <= ss <= slen or not 1 <= se <= slen
                or abs(qe - qs) + 1 > length or abs(se - ss) + 1 > length
            ):
                raise ValueError(f"Invalid BLAST score/length/coordinates: {path}:{number}")
            hit = Hit(identity, 100 * (abs(qe - qs) + 1) / qlen,
                      100 * (abs(se - ss) + 1) / slen, bitscore)
            previous = hits.get((r, p))
            rank = lambda h: (h.bitscore, min(h.rust_coverage, h.pan_coverage), h.identity,
                              h.rust_coverage, h.pan_coverage)
            if previous is None or rank(hit) > rank(previous):
                hits[r, p] = hit
    return hits


def components(
    rust: dict[str, set[Key]], pan: dict[str, set[Key]], edges: dict[tuple[str, str], set[Key]],
) -> list[tuple[list[str], list[str]]]:
    neighbors: dict[tuple[str, str], set[tuple[str, str]]] = {
        **{("rustroo", gene): set() for gene in rust},
        **{("panaroo", gene): set() for gene in pan},
    }
    for r, p in edges:
        neighbors["rustroo", r].add(("panaroo", p))
        neighbors["panaroo", p].add(("rustroo", r))
    visited = set()
    result = []
    for node in sorted(neighbors):
        if node in visited:
            continue
        visited.add(node)
        stack = [node]
        r, p = [], []
        while stack:
            current = stack.pop()
            (r if current[0] == "rustroo" else p).append(current[1])
            for neighbor in sorted(neighbors[current]):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        result.append((sorted(r), sorted(p)))
    return result


def relationship(r: list[str], p: list[str]) -> str:
    if not p:
        return "rustroo_only"
    if not r:
        return "panaroo_only"
    if len(r) == len(p) == 1:
        return "one_to_one"
    if len(p) == 1:
        return "rustroo_partition_of_panaroo"
    if len(r) == 1:
        return "rustroo_union_of_panaroo"
    return "many_to_many"


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if a or b else 0.0


def write_rows(path: Path, fields: Iterable[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def audit(rust_dir: Path, pan_dir: Path, hits_path: Optional[Path], output: Path,
          min_identity: float = 70.0, min_coverage: float = 90.0) -> dict:
    if not (0 < min_identity <= 100 and 0 < min_coverage <= 100):
        raise ValueError("Sequence identity and coverage thresholds must be in (0,100]")
    if output.exists():
        raise ValueError(f"Audit output already exists; use a fresh directory: {output}")
    paths = [rust_dir / name for name in
             ("presence.tsv", "genes.fa", "members.tsv", "member_alleles.fa")]
    paths += [pan_dir / name for name in
              ("gene_presence_absence.csv", "gene_data.csv", "pan_genome_reference.fa")]
    if hits_path is not None:
        paths.append(hits_path)
    progress("Checking input files")
    for path in paths:
        if not path.is_file():
            raise ValueError(f"Required audit input is missing or not a file: {path}")
    hashes = {}
    for path in paths:
        progress(f"Hashing {path}")
        hashes[str(path.resolve())] = digest(path)
    progress("Reading reference sequences and Rustroo presence")
    samples, rust_presence = load_presence(rust_dir / "presence.tsv")
    rust_refs = fasta_sequences(rust_dir / "genes.fa")
    if set(rust_refs) != set(rust_presence):
        raise ValueError("Rustroo references and presence families differ")
    pan_refs = fasta_sequences(pan_dir / "pan_genome_reference.fa")
    if hits_path is not None:
        progress("Validating BLAST hits before loading large membership tables")
        hits = load_hits(hits_path, rust_refs, pan_refs)
    else:
        progress("No BLAST file supplied; nonidentical reference pairs will be marked not_aligned")
        hits = {}
    progress("Reading Rustroo member alleles and memberships")
    alleles = fasta_sequences(rust_dir / "member_alleles.fa")
    members = load_members(rust_dir / "members.tsv", rust_presence, alleles)
    progress(f"Rustroo memberships complete: {len(members):,} occurrences")
    pan, pan_presence, pan_data, exceptions = load_pan(pan_dir, samples)
    pan_refound = Counter(row["gene"] for row in exceptions if row["reason"] == "refound")
    if not set(pan_refs).issubset(pan):
        raise ValueError("Panaroo reference IDs absent from membership matrix")
    progress("Comparing original annotation memberships and exact member DNA")
    rust = {gene: set() for gene in rust_presence}
    by_key: dict[Key, list[Member]] = defaultdict(list)
    loci: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    member_counts: dict[str, Counter] = defaultdict(Counter)
    for member in members:
        member_counts[member.gene][member.origin] += 1
        loci[member.gene][member.sample].add(member.locus)
        if member.annotation:
            by_key[member.sample, member.annotation].append(member)
        else:
            reason = "refound" if member.origin == "refound" else "missing_annotation_id"
            member_counts[member.gene]["missing_id"] += int(reason == "missing_annotation_id")
            exceptions.append({"tool": "rustroo", "gene": member.gene, "sample": member.sample,
                               "annotation_id": "", "reason": reason})
    pan_owners: dict[Key, set[str]] = defaultdict(set)
    for gene, keys in pan.items():
        for key in keys:
            pan_owners[key].add(gene)
    excluded = set()
    excluded_r: Counter[str] = Counter()
    excluded_p: Counter[str] = Counter()
    for key, owners in pan_owners.items():
        if len(owners) > 1:
            excluded.add(key)
            for gene in sorted(owners):
                exceptions.append({"tool": "panaroo", "gene": gene, "sample": key[0],
                                   "annotation_id": key[1], "reason": "multi_family_annotation"})
    for key, rows in by_key.items():
        if len(rows) > 1 or any(row.origin == "fusion_piece" for row in rows):
            excluded.add(key)
            for gene in sorted({row.gene for row in rows}):
                exceptions.append({"tool": "rustroo", "gene": gene, "sample": key[0],
                                   "annotation_id": key[1], "reason": "fusion_or_nonunique_annotation"})
    for gene, keys in pan.items():
        excluded_p[gene] = len(keys & excluded)
        keys.difference_update(excluded)
    for key, rows in by_key.items():
        if key not in excluded:
            rust[rows[0].gene].add(key)
        else:
            for gene in {row.gene for row in rows}:
                excluded_r[gene] += 1
    edges: dict[tuple[str, str], set[Key]] = defaultdict(set)
    mismatches, contig_mismatches = set(), set()
    for gene, keys in rust.items():
        for key in keys:
            owners = pan_owners.get(key, set())
            if not owners:
                reason = "not_in_panaroo_final_matrix" if key in pan_data else "not_in_panaroo_gene_data"
                exceptions.append({"tool": "rustroo", "gene": gene, "sample": key[0],
                                   "annotation_id": key[1], "reason": reason})
                continue
            p = next(iter(owners))
            edges[gene, p].add(key)
            rows = by_key[key]
            if any(alleles[row.allele] != pan_data[key].dna for row in rows):
                mismatches.add(key)
            if any(row.contig != pan_data[key].contig for row in rows):
                contig_mismatches.add(key)
            for reason, condition in (("source_sequence_mismatch", key in mismatches),
                                      ("source_contig_mismatch", key in contig_mismatches)):
                if condition:
                    exceptions.append({"tool": "both", "gene": gene + ";" + p, "sample": key[0],
                                       "annotation_id": key[1], "reason": reason})
    for gene, keys in pan.items():
        for key in keys:
            if key not in by_key:
                exceptions.append({"tool": "panaroo", "gene": gene, "sample": key[0],
                                   "annotation_id": key[1], "reason": "not_in_rustroo_emitted_members"})

    progress("Building family relationships and sequence reports")
    top: dict[str, tuple[float, set[str]]] = {}
    for (r, p), hit in hits.items():
        previous = top.get(r)
        if previous is None or hit.bitscore > previous[0]:
            top[r] = (hit.bitscore, {p})
        elif hit.bitscore == previous[0]:
            previous[1].add(p)
    rows: dict[str, list[dict]] = {name: [] for name in FIELDS}
    pair_support = {}
    source_disagreements = mismatches | contig_mismatches
    shared_r: dict[str, set[Key]] = defaultdict(set)
    shared_p: dict[str, set[Key]] = defaultdict(set)
    for (r, p), keys in sorted(edges.items()):
        shared_r[r].update(keys)
        shared_p[p].update(keys)
        hit = hits.get((r, p))
        if p not in pan_refs:
            status = "panaroo_reference_missing"
        elif rust_refs[r] == pan_refs[p]:
            status = "full_length_exact"
        elif hit is None:
            status = "no_reported_hit" if hits_path is not None else "not_aligned"
        else:
            status = "passes_thresholds" if hit.passes(min_identity, min_coverage) else "below_thresholds"
        exact = status == "full_length_exact"
        supported = exact or bool(hit and hit.passes(min_identity, min_coverage))
        pair_support[r, p] = supported
        rows["membership_pairs.tsv"].append({
            "rustroo_gene": r, "panaroo_gene": p, "shared_annotations": len(keys),
            "rustroo_shared_fraction": len(keys) / len(rust[r]),
            "panaroo_shared_fraction": len(keys) / len(pan[p]),
            "membership_jaccard": jaccard(rust[r], pan[p]),
            "presence_jaccard": jaccard(rust_presence[r], pan_presence[p]),
            "source_sequence_exact": len(keys - mismatches),
            "source_sequence_mismatches": len(keys & mismatches),
            "source_contig_mismatches": len(keys & contig_mismatches),
            "reference_status": status,
            "identity_pct": 100.0 if exact else hit.identity if hit else "",
            "rustroo_coverage_pct": 100.0 if exact else hit.rust_coverage if hit else "",
            "panaroo_coverage_pct": 100.0 if exact else hit.pan_coverage if hit else "",
            "passes_sequence_thresholds": int(supported),
        })
    for number, (r, p) in enumerate(components(rust, pan, edges), 1):
        component = f"c{number:06d}"
        kind = relationship(r, p)
        cooccurring = Counter(sample for gene in r for sample in rust_presence[gene])
        multicopy = {sample for gene in r for sample, ids in loci[gene].items() if len(ids) > 1}
        cooccur = sum(count > 1 for count in cooccurring.values())
        component_pairs = [(gene, target) for gene in r
                           for target in {owner for key in rust[gene] for owner in pan_owners.get(key, set())}
                           if (gene, target) in edges]
        supported_pairs = sum(pair_support[pair] for pair in component_pairs)
        if any(shared_r[gene] & source_disagreements for gene in r):
            review = "check_source_sequence_or_contig_disagreement"
        elif kind == "rustroo_partition_of_panaroo" and cooccur:
            review = "cooccurring_rustroo_groups_do_not_merge_by_count"
        elif kind == "rustroo_partition_of_panaroo" and len(component_pairs) == supported_pairs:
            review = "sequence_supported_disjoint_partition_review"
        elif kind == "rustroo_partition_of_panaroo":
            review = "partition_needs_sequence_review"
        elif kind in {"rustroo_only", "panaroo_only"}:
            review = "no_unambiguous_shared_source_annotations"
        else:
            review = "membership_and_reference_review"
        rows["components.tsv"].append({
            "component": component, "relationship": kind, "rustroo_groups": len(r),
            "panaroo_groups": len(p), "group_delta": len(r) - len(p),
            "shared_annotations": len(set().union(*(shared_r[gene] for gene in r))),
            "rustroo_cooccurring_samples": cooccur,
            "rustroo_multilocus_samples": len(multicopy),
            "sequence_supported_membership_pairs": supported_pairs,
            "membership_pairs": len(component_pairs), "review_evidence": review,
        })
        for tool, names, groups, presence, refs, shared in (
            ("rustroo", r, rust, rust_presence, rust_refs, shared_r),
            ("panaroo", p, pan, pan_presence, pan_refs, shared_p),
        ):
            for gene in names:
                keys = groups[gene]
                rows["families.tsv"].append({
                    "tool": tool, "gene": gene, "component": component, "relationship": kind,
                    "samples": len(presence[gene]), "source_annotations": len(keys),
                    "excluded_source_annotations": excluded_r[gene] if tool == "rustroo" else excluded_p[gene],
                    "shared_annotations": len(shared[gene]),
                    "shared_fraction": len(shared[gene]) / len(keys) if keys else "",
                    "unshared_annotations": len(keys - shared[gene]),
                    "refound_members": member_counts[gene]["refound"] if tool == "rustroo" else
                        pan_refound[gene],
                    "fusion_pieces": member_counts[gene]["fusion_piece"] if tool == "rustroo" else "",
                    "missing_id_members": member_counts[gene]["missing_id"] if tool == "rustroo" else "",
                    "source_sequence_mismatches": len(keys & mismatches),
                    "source_contig_mismatches": len(keys & contig_mismatches),
                    "reference_bp": len(refs[gene]) if gene in refs else "",
                    "top_sequence_hits": ";".join(sorted(top[gene][1])) if tool == "rustroo" and gene in top else "",
                    "top_hit_tied": int(len(top[gene][1]) > 1) if tool == "rustroo" and gene in top else "",
                })
    for (r, p), hit in sorted(hits.items()):
        rows["sequence_hits.tsv"].append({
            "rustroo_gene": r, "panaroo_gene": p, "identity_pct": hit.identity,
            "rustroo_coverage_pct": hit.rust_coverage, "panaroo_coverage_pct": hit.pan_coverage,
            "bitscore": hit.bitscore,
            "passes_sequence_thresholds": int(hit.passes(min_identity, min_coverage)),
            "is_top_bitscore": int(p in top[r][1]), "has_shared_annotations": int((r, p) in edges),
        })
    rows["annotation_exceptions.tsv"] = sorted(
        exceptions, key=lambda row: (row["tool"], row["gene"], row["sample"], row["annotation_id"], row["reason"]),
    )
    breakdown: dict[str, dict[str, int]] = {}
    for row in rows["components.tsv"]:
        entry = breakdown.setdefault(row["relationship"], {"components": 0, "rustroo_groups": 0,
                                                          "panaroo_groups": 0, "group_delta": 0})
        entry["components"] += 1
        for field in ("rustroo_groups", "panaroo_groups", "group_delta"):
            entry[field] += row[field]
    summary = {
        "schema_version": 1, "samples": len(samples), "rustroo_groups": len(rust),
        "panaroo_groups": len(pan), "group_delta": len(rust) - len(pan),
        "shared_source_annotations": len(set().union(*edges.values())),
        "excluded_ambiguous_source_annotations": len(excluded),
        "source_sequence_mismatches": len(mismatches),
        "source_contig_mismatches": len(contig_mismatches),
        "relationships": breakdown, "annotation_exception_counts": dict(Counter(
            row["tool"] + ":" + row["reason"] for row in exceptions)),
        "sequence_mode": "reported_blast_hits_and_exact_references" if hits_path else "exact_references_only",
        "sequence_min_identity_pct": min_identity, "sequence_min_both_coverage_pct": min_coverage,
        "membership_denominator": "unambiguous original source annotations only; excludes fusion/nonunique sources, refound occurrences and missing IDs",
        "blast_format": BLAST_FORMAT, "input_sha256": hashes,
        "limitations": [
            "Neither tool is ground truth. Family relationships describe shared annotations, not errors.",
            "Same sample-scoped IDs do not prove identical input files; source sequence/scaffold differences are reported.",
            "Fusion/multi-family sources, missing IDs and refound occurrences are excluded from membership edges.",
            "Sequence similarity never creates a membership edge. Presence Jaccard is not annotation agreement.",
            "BLAST uses the best single HSP per pair, not summed overlapping HSP coverage; tied subjects are retained.",
            "No reported hit is not proof of absent homology; BLAST search limits and thresholds matter.",
        ],
    }
    progress("Rechecking input hashes before writing reports")
    for path in paths:
        progress(f"Rechecking {path}")
        if digest(path) != hashes[str(path.resolve())]:
            raise ValueError("Audit inputs changed while being read; rerun on completed, unchanged outputs")
    progress(f"Writing reports to {output}")
    output.mkdir(parents=True)
    for name, fields in FIELDS.items():
        write_rows(output / name, fields, rows[name])
    summary["artifact_sha256"] = {name: digest(output / name) for name in FIELDS}
    summary["audit_driver_sha256"] = digest(Path(__file__))
    summary["benchmark_helpers_sha256"] = digest(Path(__file__).with_name("benchmark.py"))
    save_json(output / "summary.json", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rustroo", required=True, type=Path)
    parser.add_argument("--panaroo", required=True, type=Path)
    parser.add_argument("--hits", type=Path, help="Rustroo-query/Panaroo-target BLAST TSV, 13 columns")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-identity", type=float, default=70.0)
    parser.add_argument("--min-coverage", type=float, default=90.0)
    args = parser.parse_args()
    try:
        summary = audit(args.rustroo, args.panaroo, args.hits, args.output,
                        args.min_identity, args.min_coverage)
    except (OSError, ValueError, csv.Error) as error:
        print(f"family audit: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        progress("Interrupted; audit did not complete. If a report directory exists, use a fresh output directory.")
        return 130
    print(json.dumps({key: summary[key] for key in
                      ("rustroo_groups", "panaroo_groups", "group_delta", "relationships")}, indent=2))
    print(f"Tables and provenance: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
