"""Splitter and capture arithmetic only. No GFF, Slurm, or tool execution."""

from pathlib import Path
import tempfile
import unittest

import macs_curve


class SplitTests(unittest.TestCase):
    def test_twenty_sizes_run_from_one_to_the_stop(self) -> None:
        sizes = macs_curve.curve_sizes()
        self.assertEqual(len(sizes), 20)
        self.assertEqual(sizes[0], 1)
        self.assertEqual(sizes[-1], 1500)
        self.assertEqual(sizes, tuple(sorted(set(sizes))))
        self.assertGreater(
            sum(1 for size in sizes if size <= 500),
            sum(1 for size in sizes if size > 500),
        )

    def test_small_cohorts_are_disjoint_and_the_large_ones_reuse_genomes(self) -> None:
        names = [f"s{i:04d}.gff3" for i in range(1514)]
        first = macs_curve.split_cohorts(names, seed=1729)
        again = macs_curve.split_cohorts(names, seed=1729)
        self.assertEqual(len(first), 20)
        partials = [cohort for cohort in first if cohort["disjoint"]]
        self.assertGreaterEqual(len(partials), 5)
        seen = [name for cohort in partials for name in cohort["files"]]
        self.assertEqual(len(seen), len(set(seen)))
        self.assertEqual(sum(cohort["n"] for cohort in partials), len(seen))
        self.assertLessEqual(len(seen), len(names))
        largest = first[-1]
        self.assertFalse(largest["disjoint"])
        self.assertEqual(largest["n"], 1500)
        self.assertEqual(len(largest["files"]), 1500)
        self.assertEqual(
            [cohort["files"] for cohort in first],
            [cohort["files"] for cohort in again],
        )

    def test_one_node_keeps_the_same_cohort_as_the_full_split(self) -> None:
        names = [f"s{i:04d}.gff3" for i in range(1514)]
        cohorts = macs_curve.split_cohorts(names, seed=1729)
        one = macs_curve.choose_cohorts(cohorts, task=0)
        self.assertEqual(one[0]["n"], 1)
        self.assertEqual(one[0]["files"], cohorts[0]["files"])
        full = macs_curve.choose_cohorts(cohorts, "full")
        self.assertEqual(full[0]["n"], 1500)
        self.assertFalse(full[0]["disjoint"])
        self.assertEqual(macs_curve.choose_cohorts(cohorts, None), cohorts)

    def test_a_short_collection_caps_the_large_steps(self) -> None:
        cohorts = macs_curve.split_cohorts([f"s{i:04d}.gff3" for i in range(100)], seed=1)
        self.assertEqual(cohorts[-1]["n"], 100)
        self.assertFalse(cohorts[-1]["disjoint"])
        disjoint = [cohort for cohort in cohorts if cohort["disjoint"]]
        seen = [name for cohort in disjoint for name in cohort["files"]]
        self.assertEqual(len(seen), len(set(seen)))


class CaptureTests(unittest.TestCase):
    def test_capture_is_the_share_of_each_reference_set(self) -> None:
        fields = macs_curve.capture_fields({
            "panaroo_references": 10,
            "rustroo_references": 8,
            "accepted_pairs": 7,
            "mean_identity_pct": 99.0,
            "weighted_identity_pct": 98.5,
            "full_length_exact": 4,
        })
        self.assertEqual(fields["panaroo_captured_pct"], 70.0)
        self.assertEqual(fields["rustroo_captured_pct"], 87.5)
        self.assertEqual(fields["mean_identity_pct"], 99.0)

    def test_both_tools_receive_the_same_thread_count(self) -> None:
        gffs = [Path("a.gff3"), Path("b.gff3")]
        rust = macs_curve.rustroo_argv(Path("rustroo"), gffs, Path("out"), 32, "high")
        pan = macs_curve.panaroo_argv(Path("pan.sif"), gffs, Path("cohort"), 32, Path("gffs"))
        self.assertEqual(rust[rust.index("--threads") + 1], "32")
        self.assertEqual(rust[rust.index("--quality") + 1], "high")
        self.assertEqual(pan[pan.index("-t") + 1], "32")
        self.assertNotIn("--verbose", pan)
        self.assertNotIn("--quiet", pan)
        self.assertIn("OMP_NUM_THREADS=32", pan)

    def test_one_genome_is_a_panaroo_list_rather_than_one_gff(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cohort = Path(temporary)
            argv = macs_curve.panaroo_argv(Path("pan.sif"), [Path("sample.gff")], cohort, 4, Path("gffs"))
            self.assertEqual((cohort / "panaroo_inputs.txt").read_text(), "/inputs/sample.gff\n")
            self.assertIn("/cohort/panaroo_inputs.txt", argv)
            self.assertNotIn("/inputs/sample.gff", argv)

    def test_curve_diagram_uses_finished_points(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "n0050").mkdir()
            macs_curve.save_json(root / "n0050" / "result.json", {
                "n": 50, "disjoint": True, "status": "ok",
                "rustroo_wall_min": 1.0, "panaroo_wall_min": 4.0,
                "panaroo_captured_pct": 90.0, "rustroo_captured_pct": 88.0,
                "mean_identity_pct": 99.0,
            })
            rows = macs_curve.publish(root)
            self.assertEqual(len(rows), 1)
            diagram = (root / "curve_time.svg").read_text()
            self.assertIn(">50<", diagram)
            self.assertIn("rustroo", diagram)
            self.assertIn("panaroo", diagram)
            accuracy = (root / "curve_accuracy.svg").read_text()
            self.assertIn("mean alignment identity", accuracy)
