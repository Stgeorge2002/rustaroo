"""Synthetic tests only: no GFF processing, Docker, or real tool execution."""

import csv
from pathlib import Path
import tempfile
import unittest
import subprocess
from unittest.mock import patch

import benchmark as bench


def hit(query, subject, identity=100, span=100, score=200, reverse=False):
    start, end = (span, 1) if reverse else (1, span)
    return [
        query, subject, str(identity), str(span), str(round(identity * span / 100)),
        "100", "100", "1", str(span), str(start), str(end), "1e-30", str(score),
    ]


def write_hits(path, rows):
    with path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)


class CaptureTests(unittest.TestCase):
    def test_success_preserves_output(self):
        with patch.object(bench.subprocess, "check_output", return_value="  version 1\n"):
            self.assertEqual(bench.capture(["tool", "--version"]), "version 1")

    def test_failure_surfaces_the_captured_diagnostic(self):
        failure = subprocess.CalledProcessError(
            127, ["docker", "run"], output="bash: /usr/bin/time: No such file or directory\n",
        )
        with patch.object(bench.subprocess, "check_output", side_effect=failure):
            with self.assertRaisesRegex(RuntimeError, "/usr/bin/time: No such file"):
                bench.capture(["docker", "run"])

    def test_failure_without_output_is_still_explicit(self):
        failure = subprocess.CalledProcessError(127, ["tool"], output="")
        with patch.object(bench.subprocess, "check_output", side_effect=failure):
            with self.assertRaisesRegex(RuntimeError, "produced no diagnostic output"):
                bench.capture(["tool"])


class DesignTests(unittest.TestCase):
    def test_eight_pairs_without_warmups(self):
        names = [f"sample{i:02d}.gff3" for i in range(44)]
        sizes = list(range(5, 41, 5))
        subsets, jobs = bench.make_design(names, [1729], sizes, 1, 11729, warmups=False)
        self.assertEqual(len(subsets), 8)
        self.assertEqual(len(jobs), 16)
        self.assertTrue(all(not job["warmup"] and job["repeat"] == 1 for job in jobs))
        for count in sizes:
            selected = [job for job in jobs if job["n"] == count]
            self.assertEqual({job["tool"] for job in selected}, set(bench.TOOLS))
            self.assertEqual(len(selected), 2)
            self.assertEqual(len(subsets[f"s01_n{count:02d}"]), count)

    def test_defaults_are_nested_reproducible_and_balanced(self):
        names = [f"sample{i:02d}.gff3" for i in range(44)]
        seeds = list(range(1729, 1734))
        subsets, jobs = bench.make_design(names, seeds, list(bench.SIZES), 5, 11729)
        self.assertEqual((subsets, jobs), bench.make_design(names[::-1], seeds, list(bench.SIZES), 5, 11729))
        self.assertEqual(len(subsets), 45)
        self.assertEqual(sum(not job["warmup"] for job in jobs), 450)
        self.assertEqual(sum(job["warmup"] for job in jobs), 90)
        self.assertEqual(len({job["id"] for job in jobs}), 540)
        for series in range(1, 6):
            previous = set()
            for count in bench.SIZES:
                selected = subsets[f"s{series:02d}_n{count:02d}"]
                self.assertEqual(len(selected), count)
                self.assertEqual(selected, sorted(selected))
                self.assertTrue(previous.issubset(selected))
                previous = set(selected)
        for key in subsets:
            for tool in bench.TOOLS:
                self.assertEqual(sum(
                    job["subset"] == key and job["tool"] == tool and not job["warmup"]
                    for job in jobs
                ), 5)
        self.assertGreater(len({tuple(subsets[f"s{s:02d}_n03"]) for s in range(1, 6)}), 1)

    def test_invalid_designs(self):
        for names, seeds, sizes, repeats in (
            (["a"] * 40, [1], [3], 1), (list(map(str, range(40))), [1, 1], [3], 1),
            (list(map(str, range(40))), [], [3], 1),
            (["a", "b"], [1], [3], 1), (list(map(str, range(40))), [1], [3, 3], 1),
            (list(map(str, range(40))), [1], [3], 0),
        ):
            with self.assertRaises(ValueError):
                bench.make_design(names, seeds, sizes, repeats, 1)

    def test_timing_commands_exclude_docker_startup_and_pin_image(self):
        config = {
            "cpus": "0,1", "threads": 2, "subsets": {"s01_n03": ["a.gff3", "b.gff3", "c.gff3"]},
            "image_id": "sha256:pinned", "run_id": "test", "clean_mode": "strict",
            "tool_environment": {"RAYON_NUM_THREADS": "2", "OMP_NUM_THREADS": "2"},
        }
        job = {"id": "testjob", "subset": "s01_n03", "tool": "panaroo"}
        with patch.object(bench.os, "getuid", return_value=1000, create=True), \
             patch.object(bench.os, "getgid", return_value=1000, create=True):
            pan = bench.tool_command(config, Path("/benchmark"), job)
        self.assertLess(pan.index("sha256:pinned"), pan.index("/usr/bin/time"))
        self.assertIn("--pull=never", pan)
        self.assertIn("--cpuset-cpus", pan)
        self.assertIn("--remove-invalid-genes", pan)
        self.assertIn("RAYON_NUM_THREADS=2", pan)
        native = bench.tool_command(config, Path("/benchmark"), {**job, "tool": "rustroo"})
        self.assertEqual(native[:3], ["taskset", "--cpu-list", "0,1"])
        self.assertEqual(native[native.index("--threads") + 1], "2")


class MetricsTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)

    def tearDown(self):
        self.directory.cleanup()

    def test_timing_metrics_and_bad_values(self):
        path = self.root / "timing.txt"
        valid = "wall_seconds=10\nuser_seconds=20\nsystem_seconds=5\nmax_process_rss_kib=1024\nexit_code=0\n"
        path.write_text(valid)
        self.assertEqual(bench.timing_metrics(path)["average_cpu_cores"], 2.5)
        for text in (valid.replace("exit_code=0", "exit_code=1"), valid.replace("=10", "=0"),
                     valid.replace("=20", "=nan"), "wall_seconds=10\n"):
            path.write_text(text)
            with self.assertRaises(ValueError):
                bench.timing_metrics(path)

    def test_fasta_and_presence_counts_and_sample_validation(self):
        output = self.root / "output"
        output.mkdir()
        (self.root / "tool.log").write_text("invalid gene! file - id: sample - cds\n")
        (output / "genes.fa").write_text(">a\nATG\n>b\nAAA\n")
        (output / "presence.tsv").write_text("gene\ts1\ts2\ts3\na\t1\t1\t1\nb\t1\t1\t0\n")
        names = ["s1.gff3", "s2.gff3", "s3.gff3"]
        result = bench.output_metrics(output, "rustroo", names)
        self.assertEqual(result["groups"], 2)
        self.assertEqual(result["references"], 2)
        self.assertEqual(result["core_groups_99pct"], 1)
        self.assertEqual(result["invalid_cds_log_messages"], 1)
        (output / "genes.fa").write_text(">b\nAAA\n>a\nATG\n")
        self.assertEqual(bench.output_metrics(output, "rustroo", names)["reference_sha256"], result["reference_sha256"])
        with self.assertRaises(ValueError):
            bench.output_metrics(output, "rustroo", ["other.gff3"] * 3)
        (output / "genes.fa").write_text(">a\nATG\n>a\nAAA\n")
        with self.assertRaises(ValueError):
            bench.fasta_sequences(output / "genes.fa")

    def test_reciprocity_ties_coverage_reverse_strand_and_empty_matches(self):
        forward, reverse = self.root / "forward.tsv", self.root / "reverse.tsv"
        write_hits(forward, [
            hit("A", "a", reverse=True), hit("B", "b", span=89),
            hit("C", "c", identity=99, span=95),
            hit("D", "d1"), hit("D", "d2"),
        ])
        write_hits(reverse, [hit("a", "A"), hit("b", "B"), hit("c", "C"), hit("d1", "D")])
        pan = {identifier: "A" * 100 for identifier in ("A", "B", "C", "D", "E")}
        rust = {identifier: "A" * 100 for identifier in ("a", "b", "c", "d1", "d2", "e")}
        metrics, pairs = bench.compare_references(forward, reverse, pan, rust, 70, 90)
        self.assertEqual(metrics["accepted_pairs"], 2)
        self.assertEqual(metrics["reciprocal_pairs"], 3)
        self.assertEqual(metrics["reciprocal_filtered_out"], 1)
        self.assertEqual(metrics["full_length_exact"], 1)
        self.assertEqual(metrics["panaroo_ambiguous"], 1)
        self.assertEqual(metrics["panaroo_no_hit"], 1)
        self.assertEqual(metrics["mean_identity_pct"], 99.5)
        self.assertEqual(len(pairs), 2)
        for tool, total in (("panaroo", len(pan)), ("rustroo", len(rust))):
            self.assertEqual(
                metrics["accepted_pairs"] + metrics["reciprocal_filtered_out"]
                + metrics[f"{tool}_no_hit"] + metrics[f"{tool}_ambiguous"]
                + metrics[f"{tool}_unique_nonreciprocal"], total,
            )
        forward.write_text("")
        reverse.write_text("")
        metrics, pairs = bench.compare_references(forward, reverse, pan, rust, 70, 90)
        self.assertIsNone(metrics["mean_identity_pct"])
        self.assertEqual(pairs, [])

    def test_tied_hsps_for_same_subject_are_not_ambiguous(self):
        path = self.root / "hits.tsv"
        write_hits(path, [hit("A", "a", span=95), hit("A", "a", span=100)])
        unique, ambiguous, _ = bench.best_hits(path)
        self.assertEqual(unique["A"][8], "100")
        self.assertFalse(ambiguous)

    def test_malformed_blast_is_explicit_error(self):
        path = self.root / "hits.tsv"
        path.write_text("not\tblast\n")
        with self.assertRaises(ValueError):
            bench.best_hits(path)
        row = hit("A", "a")
        row[7] = "0"
        write_hits(path, [row])
        with self.assertRaises(ValueError):
            bench.best_hits(path)

    def test_report_excludes_warmups_and_aggregates_subset_medians(self):
        names = [f"s{i}.gff3" for i in range(5)]
        subsets, jobs = bench.make_design(names, [1, 2], [3, 5], 2, 1)
        config = {"subsets": subsets, "jobs": jobs, "sizes": [3, 5], "seeds": [1, 2], "repeats": 2}
        for job in jobs:
            folder = self.root / "jobs" / job["id"]
            folder.mkdir(parents=True)
            wall = (10 if job["tool"] == "rustroo" else 20) * job["series"]
            bench.save_json(folder / "result.json", {
                "status": "ok", "started_utc": "2026-01-01T00:00:00Z",
                "launcher_wall_seconds": wall + 1,
                "timing": {"wall_seconds": wall, "max_process_rss_kib": 1024, "average_cpu_cores": 2},
                "outputs": {"groups": 2, "references": 2, "reference_sha256": "same", "presence_file_sha256": "same"},
            })
        for subset in subsets:
            folder = self.root / "concordance" / subset
            folder.mkdir(parents=True)
            bench.save_json(folder / "result.json", {
                "status": "ok", "metrics": {
                    "panaroo_accepted_fraction": 1, "rustroo_accepted_fraction": 1,
                    "accepted_pairs": 2, "full_length_exact": 2,
                    "mean_identity_pct": 100, "weighted_identity_pct": 100,
                },
            })
        bench.report(self.root, config)
        with (self.root / "report" / "runs.tsv").open() as handle:
            raw = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(raw), 16)
        self.assertTrue(all(row["warmup"] == "False" for row in raw))
        with (self.root / "report" / "scaling_summary.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        native = next(row for row in rows if row["tool"] == "rustroo")
        self.assertEqual(float(native["wall_seconds_median"]), 15)
        with (self.root / "report" / "speedup_summary.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(float(rows[0]["median"]), 2)
        self.assertIn("<svg", (self.root / "report" / "runtime.svg").read_text())
        (self.root / "jobs" / jobs[0]["id"] / "result.json").unlink()
        with self.assertRaises(FileNotFoundError):
            bench.report(self.root, config)

    def test_run_failure_is_retained_not_converted_to_success(self):
        job = {
            "id": "s01_n03_r01_rustroo", "subset": "s01_n03", "series": 1,
            "seed": 1, "n": 3, "repeat": 1, "tool": "rustroo", "warmup": False,
        }
        config = {
            "subsets": {"s01_n03": ["s1.gff3", "s2.gff3", "s3.gff3"]},
            "cpus": "0", "threads": 1, "tool_environment": {}, "run_id": "test",
        }
        (self.root / "jobs").mkdir()
        proc = self.root / "loadavg"
        proc.write_text("0.1 0.1 0.1 1/10 42")
        original_read = Path.read_text

        def read_text(path, *args, **kwargs):
            if str(path) == "/proc/loadavg":
                return original_read(proc, *args, **kwargs)
            return original_read(path, *args, **kwargs)

        with patch.object(bench.Path, "read_text", read_text), \
             patch.object(bench, "execute", side_effect=RuntimeError("synthetic failure")):
            with self.assertRaisesRegex(RuntimeError, "synthetic failure"):
                bench.run_job(self.root, config, job)
        result = bench.load_json(self.root / "jobs" / job["id"] / "result.json")
        self.assertEqual(result["status"], "failed")
        self.assertNotIn("timing", result)
        with self.assertRaisesRegex(ValueError, "Failed attempt"):
            bench.run_job(self.root, config, job)

    def test_eight_pair_report_needs_no_concordance_artifacts(self):
        sizes = list(range(5, 41, 5))
        names = [f"s{i}.gff3" for i in range(44)]
        subsets, jobs = bench.make_design(names, [1729], sizes, 1, 11729, warmups=False)
        config = {
            "subsets": subsets, "jobs": jobs, "sizes": sizes, "seeds": [1729],
            "repeats": 1, "concordance": False,
        }
        for job in jobs:
            folder = self.root / "jobs" / job["id"]
            folder.mkdir(parents=True)
            wall = job["n"] * (1 if job["tool"] == "rustroo" else 2)
            bench.save_json(folder / "result.json", {
                "status": "ok", "started_utc": "2026-01-01T00:00:00Z",
                "launcher_wall_seconds": wall + 1,
                "timing": {"wall_seconds": wall, "max_process_rss_kib": 1024, "average_cpu_cores": 2},
                "outputs": {
                    "groups": 2, "references": 2, "reference_sha256": "same",
                    "presence_file_sha256": "same",
                },
            })
        bench.report(self.root, config)
        with (self.root / "report" / "runs.tsv").open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 16)
        with (self.root / "report" / "paired_speedups.tsv").open() as handle:
            pairs = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(pairs), 8)
        self.assertEqual({int(row["n"]) for row in pairs}, set(sizes))
        self.assertTrue(all(float(row["panaroo_over_rustroo"]) == 2 for row in pairs))
        self.assertFalse((self.root / "concordance").exists())
        self.assertFalse((self.root / "report" / "concordance.tsv").exists())
        self.assertFalse((self.root / "report" / "matched_references.svg").exists())
        figure = (self.root / "report" / "runtime.svg").read_text()
        self.assertNotIn('text-anchor="middle">3</text>', figure)
        self.assertIn('text-anchor="middle">5</text>', figure)


if __name__ == "__main__":
    unittest.main()
