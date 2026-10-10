"""Synthetic audit tests; no Docker, BLAST or real dataset execution."""

from __future__ import annotations

import csv
from contextlib import redirect_stderr
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from audit_families import FIELDS, audit, load_hits, load_pan, main, write_rows


class FamilyAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.rust = self.root / "rust"
        self.pan = self.root / "pan"
        self.rust.mkdir()
        self.pan.mkdir()
        self.samples = ["s1", "s2"]
        self.dna = "ATGAAATAA"

    def fixture(
        self, rust: dict[str, list[tuple[str, str]]],
        pan: dict[str, list[tuple[str, str]]],
        origins: dict[tuple[str, str], str] | None = None,
    ) -> None:
        origins = origins or {}
        with (self.rust / "presence.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["gene", *self.samples])
            for gene, members in rust.items():
                present = {sample for sample, _ in members}
                writer.writerow([gene, *[int(sample in present) for sample in self.samples]])
        (self.rust / "genes.fa").write_text(
            "".join(f">{gene}\n{self.dna}\n" for gene in rust),
        )
        (self.rust / "member_alleles.fa").write_text(f">dna_0\n{self.dna}\n")
        fields = ("gene", "sample", "contig", "annotation_id", "start", "end", "strand",
                  "locus", "origin", "allele")
        rows = []
        index = 0
        for gene, members in rust.items():
            for sample, annotation in members:
                index += 1
                rows.append(dict(zip(fields, (
                    gene, sample, "ctg", annotation, index * 100, index * 100 + 8, 1,
                    str(index), origins.get((sample, annotation), "annotated"), "dna_0",
                ))))
        write_rows(self.rust / "members.tsv", fields, rows)
        with (self.pan / "gene_presence_absence.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["Gene", "Non-unique Gene name", "Annotation", *self.samples])
            for gene, members in pan.items():
                writer.writerow([gene, "", "", *[
                    ";".join(annotation for sample, annotation in members if sample == name)
                    for name in self.samples
                ]])
        with (self.pan / "gene_data.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["gff_file", "annotation_id", "scaffold_name", "dna_sequence"])
            keys = {key for members in pan.values() for key in members if "_refound_" not in key[1]}
            for sample, annotation in sorted(keys):
                writer.writerow([sample, annotation, "ctg", self.dna])
        (self.pan / "pan_genome_reference.fa").write_text(
            "".join(f">{gene}\n{self.dna}\n" for gene in pan),
        )

    def run_audit(self, hits: Path | None = None) -> dict:
        return audit(self.rust, self.pan, hits, self.root / "audit")

    def rows(self, name: str) -> list[dict[str, str]]:
        with (self.root / "audit" / name).open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_annotation_ids_are_sample_scoped(self) -> None:
        self.fixture({"r1": [("s1", "id")], "r2": [("s2", "id")]},
                     {"p1": [("s1", "id")], "p2": [("s2", "id")]})
        summary = self.run_audit()
        self.assertEqual(summary["relationships"]["one_to_one"]["components"], 2)
        self.assertEqual(summary["shared_source_annotations"], 2)

    def test_partition_accounts_for_excess_and_reports_disjoint_sequence_support(self) -> None:
        self.fixture({"r1": [("s1", "a")], "r2": [("s2", "b")]},
                     {"p": [("s1", "a"), ("s2", "b")]})
        summary = self.run_audit()
        self.assertEqual(summary["group_delta"], 1)
        component = self.rows("components.tsv")[0]
        self.assertEqual(component["relationship"], "rustroo_partition_of_panaroo")
        self.assertEqual(component["review_evidence"], "sequence_supported_disjoint_partition_review")
        pairs = self.rows("membership_pairs.tsv")
        self.assertEqual({row["reference_status"] for row in pairs}, {"full_length_exact"})
        self.assertEqual({row["membership_jaccard"] for row in pairs}, {"0.5"})

    def test_cooccurring_partition_is_not_a_merge_recommendation(self) -> None:
        self.fixture({"r1": [("s1", "a")], "r2": [("s1", "b")]},
                     {"p": [("s1", "a"), ("s1", "b")]})
        self.run_audit()
        row = self.rows("components.tsv")[0]
        self.assertEqual(row["rustroo_cooccurring_samples"], "1")
        self.assertEqual(row["review_evidence"], "cooccurring_rustroo_groups_do_not_merge_by_count")

    def test_union_many_to_many_and_orphans_decompose_total_delta(self) -> None:
        self.fixture(
            {"r1": [("s1", "a"), ("s1", "b")], "r2": [("s2", "c"), ("s2", "d")],
             "r3": [("s1", "x")], "r4": [("s2", "y")], "r5": [("s1", "only_r")]},
            {"p1": [("s1", "a"), ("s2", "c")], "p2": [("s1", "b"), ("s2", "d")],
             "p3": [("s1", "x")], "p4": [("s2", "y")], "p5": [("s2", "only_p")]},
        )
        # r3 now contains members of two Panaroo groups.
        matrix = self.rust / "members.tsv"
        matrix.write_text(matrix.read_text().replace("r4\ts2", "r3\ts2"))
        (self.rust / "presence.tsv").write_text("gene\ts1\ts2\nr1\t1\t0\nr2\t0\t1\nr3\t1\t1\nr5\t1\t0\n")
        (self.rust / "genes.fa").write_text(
            "".join(f">{gene}\n{self.dna}\n" for gene in ["r1", "r2", "r3", "r5"]),
        )
        summary = self.run_audit()
        kinds = summary["relationships"]
        self.assertIn("many_to_many", kinds)
        self.assertIn("rustroo_union_of_panaroo", kinds)
        self.assertIn("rustroo_only", kinds)
        self.assertIn("panaroo_only", kinds)
        self.assertEqual(sum(row["group_delta"] for row in kinds.values()), summary["group_delta"])

    def test_fusion_and_duplicate_family_sources_do_not_create_false_edges(self) -> None:
        self.fixture({"r1": [("s1", "fusion")], "r2": [("s1", "fusion")]},
                     {"p": [("s1", "fusion")]}, {("s1", "fusion"): "fusion_piece"})
        summary = self.run_audit()
        self.assertEqual(summary["shared_source_annotations"], 0)
        self.assertEqual(summary["excluded_ambiguous_source_annotations"], 1)
        self.assertEqual(self.rows("membership_pairs.tsv"), [])
        self.assertTrue(all(row["excluded_source_annotations"] == "1" for row in self.rows("families.tsv")))
        self.assertEqual(sum(row["group_delta"] for row in summary["relationships"].values()), 1)

    def test_refound_and_missing_ids_are_explicit_exceptions(self) -> None:
        self.fixture({"r": [("s1", "")]}, {"p": [("s1", "0_refound_1")]},
                     {("s1", ""): "refound"})
        summary = self.run_audit()
        self.assertEqual(summary["shared_source_annotations"], 0)
        self.assertEqual(summary["annotation_exception_counts"]["rustroo:refound"], 1)
        self.assertEqual(summary["annotation_exception_counts"]["panaroo:refound"], 1)
        self.assertEqual(len(self.rows("annotation_exceptions.tsv")), 2)

    def test_identical_references_do_not_invent_membership_edges(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "b")]})
        summary = self.run_audit()
        self.assertEqual(summary["shared_source_annotations"], 0)
        self.assertEqual(self.rows("membership_pairs.tsv"), [])

    def test_source_sequence_mismatch_is_visible_even_when_references_are_exact(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        data = self.pan / "gene_data.csv"
        data.write_text(data.read_text().replace(self.dna, "ATGCCCTAA"))
        summary = self.run_audit()
        self.assertEqual(summary["source_sequence_mismatches"], 1)
        pair = self.rows("membership_pairs.tsv")[0]
        self.assertEqual(pair["source_sequence_exact"], "0")
        self.assertEqual(pair["reference_status"], "full_length_exact")
        self.assertEqual(self.rows("components.tsv")[0]["review_evidence"],
                         "check_source_sequence_or_contig_disagreement")

    def test_source_contig_mismatch_is_visible(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        data = self.pan / "gene_data.csv"
        data.write_text(data.read_text().replace(",ctg,", ",different,"))
        summary = self.run_audit()
        self.assertEqual(summary["source_contig_mismatches"], 1)

    def test_panaroo_dna_pool_preserves_exact_case_normalized_sequences(self) -> None:
        self.fixture({"r": [("s1", "a"), ("s2", "b")]},
                     {"p": [("s1", "a"), ("s2", "b")]})
        path = self.pan / "gene_data.csv"
        with path.open("a", newline="") as handle:
            csv.writer(handle).writerow(["s1", "c", "ctg", "ATGCCCTAA"])
        path.write_text(path.read_text().replace(self.dna, self.dna.lower()))
        _, _, data, _ = load_pan(self.pan, self.samples)
        self.assertEqual(data["s1", "a"].dna, self.dna)
        self.assertIs(data["s1", "a"].dna, data["s2", "b"].dna)
        self.assertIs(data["s1", "a"].contig, data["s2", "b"].contig)
        self.assertNotEqual(data["s1", "a"].dna, data["s1", "c"].dna)

    def test_missing_hits_fail_before_loading_large_panaroo_table(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        with patch("audit_families.load_pan") as reader:
            with self.assertRaisesRegex(ValueError, "Required audit input"):
                self.run_audit(self.root / "missing_hits.tsv")
            reader.assert_not_called()
        self.assertFalse((self.root / "audit").exists())

    def test_malformed_hits_fail_before_loading_large_panaroo_table(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        hits = self.root / "hits.tsv"
        hits.write_text("malformed\n")
        with patch("audit_families.load_pan") as reader:
            with self.assertRaisesRegex(ValueError, "Expected"):
                self.run_audit(hits)
            reader.assert_not_called()
        self.assertFalse((self.root / "audit").exists())

    def test_progress_makes_slow_stages_and_hash_recheck_visible(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        messages = io.StringIO()
        with redirect_stderr(messages):
            self.run_audit()
        text = messages.getvalue()
        for stage in ("Checking input files", "Hashing", "No BLAST file supplied",
                      "Reading Panaroo gene_data.csv", "1 annotations, 1 distinct DNA sequences",
                      "Comparing original annotation", "Rechecking input hashes", "Writing reports"):
            self.assertIn(stage, text)

    def test_interrupt_returns_130_and_explains_incomplete_audit(self) -> None:
        messages = io.StringIO()
        with patch("sys.argv", ["audit_families.py", "--rustroo", str(self.rust),
                                "--panaroo", str(self.pan), "--output", str(self.root / "audit")]):
            with patch("audit_families.audit", side_effect=KeyboardInterrupt), redirect_stderr(messages):
                self.assertEqual(main(), 130)
        self.assertIn("audit did not complete", messages.getvalue())

    def test_unknown_panaroo_ids_fail_instead_of_being_guessed(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        (self.pan / "gene_data.csv").write_text("gff_file,annotation_id,scaffold_name,dna_sequence\n")
        with self.assertRaisesRegex(ValueError, "absent from gene_data"):
            self.run_audit()

    def test_mismatched_cohorts_fail(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        path = self.pan / "gene_presence_absence.csv"
        path.write_text(path.read_text().replace("s2", "wrong_sample"))
        with self.assertRaisesRegex(ValueError, "sample sets differ"):
            self.run_audit()

    def test_presence_member_inconsistency_fails(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        (self.rust / "presence.tsv").write_text("gene\ts1\ts2\nr\t1\t1\n")
        with self.assertRaisesRegex(ValueError, "sample sets differ"):
            self.run_audit()

    def test_missing_panaroo_reference_is_not_a_missing_family(self) -> None:
        self.fixture({"r1": [("s1", "a")], "r2": [("s2", "b")]},
                     {"p1": [("s1", "a")], "p2": [("s2", "b")]})
        (self.pan / "pan_genome_reference.fa").write_text(f">p1\n{self.dna}\n")
        summary = self.run_audit()
        self.assertEqual(summary["panaroo_groups"], 2)
        pairs = self.rows("membership_pairs.tsv")
        self.assertEqual(pairs[1]["reference_status"], "panaroo_reference_missing")

    def test_tied_sequence_hits_are_retained_without_creating_membership(self) -> None:
        self.fixture({"r": [("s1", "a")]},
                     {"p1": [("s1", "a")], "p2": [("s2", "b")]})
        hits = self.root / "hits.tsv"
        hits.write_text(
            "r\tp1\t100\t9\t9\t9\t9\t1\t9\t1\t9\t1e-10\t20\n"
            "r\tp2\t100\t9\t9\t9\t9\t1\t9\t9\t1\t1e-10\t20\n",
        )
        self.run_audit(hits)
        sequences = self.rows("sequence_hits.tsv")
        self.assertEqual(len(sequences), 2)
        self.assertEqual([row["is_top_bitscore"] for row in sequences], ["1", "1"])
        self.assertEqual([row["has_shared_annotations"] for row in sequences], ["1", "0"])
        family = next(row for row in self.rows("families.tsv") if row["tool"] == "rustroo")
        self.assertEqual(family["top_hit_tied"], "1")
        self.assertEqual(family["top_sequence_hits"], "p1;p2")

    def test_single_hsp_coverage_is_not_summed_and_boundary_is_inclusive(self) -> None:
        path = self.root / "hits.tsv"
        path.write_text(
            "r\tp\t70\t10\t7\t10\t10\t1\t9\t1\t9\t1e-10\t20\n"
            "r\tp\t100\t5\t5\t10\t10\t6\t10\t6\t10\t1e-10\t10\n",
        )
        hit = load_hits(path, {"r": "A" * 10}, {"p": "A" * 10})["r", "p"]
        self.assertEqual(hit.rust_coverage, 90.0)
        self.assertTrue(hit.passes(70, 90))
        self.assertFalse(hit.passes(70, 90.01))

    def test_invalid_blast_is_rejected(self) -> None:
        path = self.root / "hits.tsv"
        for row in (
            "r\tp\t100\t9\t9\t9\t9\t1\t10\t1\t9\t1e-10\t20",
            "r\tp\t100\t9\t9\t9\t9\t1\t9\t1\t9\t1e-10\tnan",
            "r\tp\t100\t9\t9\t9\t9\t1\t9\t1\t9\t-1\t20",
            "wrong\tp\t100\t9\t9\t9\t9\t1\t9\t1\t9\t1e-10\t20",
            "malformed",
        ):
            path.write_text(row + "\n")
            with self.assertRaises(ValueError):
                load_hits(path, {"r": self.dna}, {"p": self.dna})

    def test_empty_tables_have_headers_and_completed_audit_cannot_be_overwritten(self) -> None:
        self.fixture({"r": [("s1", "a")]}, {"p": [("s1", "a")]})
        summary = self.run_audit()
        for name, fields in FIELDS.items():
            self.assertTrue((self.root / "audit" / name).read_text().startswith("\t".join(fields)))
        self.assertEqual(len(summary["input_sha256"]), 7)
        self.assertEqual(len(summary["artifact_sha256"]), len(FIELDS))
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.run_audit()


if __name__ == "__main__":
    unittest.main()
