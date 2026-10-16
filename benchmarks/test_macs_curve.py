"""Splitter and capture arithmetic only. No GFF, Slurm, or tool execution."""

from pathlib import Path
import tempfile
import unittest

import macs_curve


class SplitTests(unittest.TestCase):
    def test_partial_cohorts_are_disjoint_and_the_full_set_is_everyone(self) -> None:
        names = [f"s{i:04d}.gff3" for i in range(1400)]
        first = macs_curve.split_cohorts(names, seed=1729)
        again = macs_curve.split_cohorts(names, seed=1729)
        partials = [cohort for cohort in first if cohort["disjoint"]]
        self.assertEqual([cohort["n"] for cohort in partials], [50, 150, 300, 800])
        seen = [name for cohort in partials for name in cohort["files"]]
        self.assertEqual(len(seen), len(set(seen)))
        self.assertEqual(sum(cohort["n"] for cohort in partials), len(seen))
        full = first[-1]
        self.assertFalse(full["disjoint"])
        self.assertEqual(full["n"], 1400)
        self.assertEqual(set(full["files"]), set(names))
        self.assertEqual(
            [cohort["files"] for cohort in first],
            [cohort["files"] for cohort in again],
        )

    def test_too_few_assemblies_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            macs_curve.split_cohorts(["a.gff3"] * 100)


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
        rust = macs_curve.rustroo_argv(Path("rustroo"), gffs, Path("out"), 32)
        pan = macs_curve.panaroo_argv(Path("pan.sif"), gffs, Path("cohort"), 32, Path("gffs"))
        self.assertEqual(rust[rust.index("--threads") + 1], "32")
        self.assertEqual(pan[pan.index("-t") + 1], "32")
        self.assertIn("OMP_NUM_THREADS=32", pan)

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
