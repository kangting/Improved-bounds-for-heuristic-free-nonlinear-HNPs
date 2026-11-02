"""Tests of JSON aggregation only: no imports from the attack modules."""
import csv
import io
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_results as audit


class AuditTests(unittest.TestCase):
    def row(self, status, elapsed, singleton=None, nodes=None):
        return {"family": "mihnp", "method": "ours", "observation_model": "independent",
            "unknown_bits": 3, "N": 3, "source_observations_used": 3,
            "source_rank": 6, "ambient_dimension": 6,
            "audit.revealed_fraction": .5, "status": status,
            "audit.algorithm_finished": status in audit.FINISHED,
            "audit.singleton_true_recounted": singleton,
            "audit.root_node_completed": nodes == 1 if status == "COMPLETE" and nodes is not None else None,
            "parent_elapsed_seconds": elapsed, "nodes_visited": nodes}

    def test_censored_times_not_mixed_with_finished(self):
        row = audit.summarize([self.row("COMPLETE", 2, True, 1),
                               self.row("CENSORED_WORK", 100, None, 4096)])[0]
        self.assertEqual(row["finished_parent_elapsed_seconds_median"], 2)
        self.assertEqual(row["interrupted_parent_elapsed_seconds_max"], 100)
        self.assertEqual(row["CENSORED_WORK_count"], 1)

    def test_no_candidate_is_finished_but_not_success(self):
        row = audit.summarize([self.row("NO_CANDIDATE", 4, False)])[0]
        self.assertEqual((row["algorithm_finished"], row["other_finished"], row["interrupted"]), (1, 1, 0))

    def test_root_not_assumed_grid_cell(self):
        row = audit.summarize([self.row("COMPLETE", 2, True, 1)])[0]
        self.assertEqual(row["root_node_completed_count"], 1)
        self.assertIsNone(row["first_grid_cell_completed_count"])
        self.assertIsNone(row["distinct_tested_targets_recorded_median"])

    def test_unrecorded_counter_is_not_zero(self):
        row = audit.summarize([self.row("NO_CANDIDATE", 3, False)])[0]
        self.assertIsNone(row["finished_decoder_calls_median"])
        self.assertIsNone(row["root_node_completed_count"])

    def test_null_output_not_empty_list(self):
        self.assertEqual(audit.recount(None, 7), (None, None, None))
        self.assertEqual(audit.recount([], 7), (0, False, False))

    def test_missing_zero_false_serialization(self):
        raw = audit.csv_bytes([{"missing": None, "zero": 0, "false": False, "empty": []}])
        row = next(csv.DictReader(io.StringIO(raw.decode())))
        self.assertEqual(row, {"empty": "[]", "false": "false", "missing": "NA", "zero": "0"})

    def test_duplicate_ids_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            audit.keyed([{"id": "a"}, {"id": "a"}], "fixture")

    def test_existing_study_read_only(self):
        files = audit.artifacts(audit.ROOT)
        summary = list(csv.DictReader(io.StringIO(files["summary.csv"].decode())))
        self.assertEqual(len(summary), 35)
        self.assertEqual(sum(int(row["runs"]) for row in summary), 1050)
        self.assertEqual(sum(int(row["singleton_true"]) for row in summary), 547)
        self.assertEqual(sum(int(row["algorithm_finished"]) for row in summary), 900)
        self.assertEqual(sum(int(row["CENSORED_WORK_count"]) for row in summary), 150)
        self.assertTrue(all(int(row["runs"]) == 30 for row in summary))


if __name__ == "__main__":
    unittest.main()
