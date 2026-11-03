"""Check stored data and verification logic; do not execute archival solvers."""
import copy
import gzip
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import verify_public as public
from evaluate_saved import evaluate


class SupplementChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reports = public.check_public(ROOT)

    def model(self, identifier=public.IDS[1]):
        job = ROOT / "records/trials" / identifier
        with gzip.open(job / "certificate.json.gz", "rt") as handle:
            certificate = json.load(handle)
        return public.load(job / "public.json"), certificate

    def test_all_three_original_certificates(self):
        self.assertEqual(set(self.reports), set(public.IDS))
        self.assertEqual([self.reports[i]["verification_replay"]["visited_nodes"] for i in public.IDS], [1, 5, 1])

    def test_basis_dimensions_and_full_identity(self):
        for report in self.reports.values():
            b = report["basis_evidence"]
            self.assertEqual((len(b["B0_reconstructed_now"]), len(b["B0_reconstructed_now"][0])), (11, 16))
            self.assertEqual(abs(b["det_U"]), 1)
            self.assertTrue(report["stored_outputs_publicly_validated"])

    def test_exact_ambiguous_root_target(self):
        root = self.reports[public.IDS[1]]["verification_replay"]["nodes"][0]
        self.assertEqual(root["final_target_interval"],
            [64107237951601111634317879330044870, 239763353792200743855781330110940573])
        self.assertEqual(root["final_target_cardinality"], 175656115840599632221463450780895704)
        self.assertEqual(root["ambiguous_quotient_rows"], [9])

    def test_five_node_split_and_pruning_trace(self):
        replay = self.reports[public.IDS[1]]["verification_replay"]
        self.assertEqual([row["decision"] for row in replay["nodes"]],
                         ["SPLIT", "EMPTY_TARGET_INTERVAL", "SPLIT", "SINGLETON_TARGET", "EMPTY_TARGET_INTERVAL"])
        self.assertEqual([row["split_coordinate0"] for row in replay["nodes"] if row["decision"] == "SPLIT"], [0, 0])
        self.assertEqual(replay["empty_target_interval_subtrees"], 2)
        self.assertEqual(replay["pending_frontier_nodes"], 0)
        self.assertEqual(replay["logical_cells_covered"], 17**11)

    def test_no_decoder_execution_or_new_times_in_replay(self):
        for report in self.reports.values():
            replay = report["verification_replay"]
            self.assertFalse(replay["decoder_executed"])
            self.assertFalse(replay["exhaustive_target_scan_executed"])
            self.assertIsNone(replay["new_timing_measurements"])
            self.assertEqual(len(replay["distinct_singleton_targets"]), 1)

    def test_nonunimodular_transform_rejected(self):
        data, certificate = self.model()
        certificate = copy.deepcopy(certificate)
        certificate["transform"][0] = [2 * value for value in certificate["transform"][0]]
        with self.assertRaisesRegex(public.VerificationError, "unimodular"):
            public.build_model(data, certificate)

    def test_wrong_target_identity_rejected(self):
        data, certificate = self.model()
        certificate = copy.deepcopy(certificate)
        certificate["alpha"][0] += 1
        with self.assertRaisesRegex(public.VerificationError, "alpha U"):
            public.build_model(data, certificate)

    def test_signed_floor_sum_without_large_enumeration(self):
        for n in (0, 1, 5):
            for a in (-3, 0, 4):
                self.assertEqual(public.floor_sum(n, 7, a, -9), sum((a * i - 9) // 7 for i in range(n)))

    def test_public_checker_never_opens_evaluation_data(self):
        original = Path.open
        evaluation = (ROOT / "evaluation").resolve()
        def guarded(path, *args, **kwargs):
            if path.resolve().is_relative_to(evaluation):
                raise AssertionError("Public checker read evaluation/truth")
            return original(path, *args, **kwargs)
        with mock.patch.object(Path, "open", guarded):
            reports = public.check_public(ROOT)
        self.assertEqual(reports, self.reports)

    def test_saved_output_constraint_check_rejects_wrong_candidate(self):
        data, _ = self.model()
        self.assertFalse(public.valid_stored_output(data, 0))
        self.assertFalse(public.valid_stored_output(data, data["reports"][0] - 1))

    def test_payload_hash_tampering_rejected(self):
        with tempfile.TemporaryDirectory(prefix="hnp_supplement_hash_") as temporary:
            folder = Path(temporary) / "copy"
            shutil.copytree(ROOT, folder, ignore=shutil.ignore_patterns("__pycache__"))
            path = folder / "records/trials" / public.IDS[0] / "public.json"
            path.write_bytes(path.read_bytes() + b" ")
            with self.assertRaisesRegex(public.VerificationError, "hash mismatch"):
                public.verify_integrity(folder)

    def test_published_trace_tampering_rejected_by_recomputation(self):
        with tempfile.TemporaryDirectory(prefix="hnp_supplement_trace_") as temporary:
            folder = Path(temporary) / "copy"
            shutil.copytree(ROOT, folder, ignore=shutil.ignore_patterns("__pycache__"))
            path = folder / "verification_replay" / (public.IDS[1] + ".json")
            value = public.load(path)
            value["verification_replay"]["nodes"][0]["decision"] = "SINGLETON_TARGET"
            path.write_bytes(public.json_bytes(value))
            with self.assertRaisesRegex(public.VerificationError, "replay differs"):
                public.verify_trial(folder, public.IDS[1])  # mathematical check, not just a hash test

    def test_replay_budget_never_expands_to_new_attack(self):
        data, certificate = self.model()
        model = public.build_model(data, certificate)
        with self.assertRaisesRegex(public.VerificationError, "stored node count"):
            public.replay(model, 1)

    def test_truth_is_a_separate_post_public_evaluation(self):
        checked = evaluate(ROOT, self.reports)
        self.assertEqual(len(checked), 3)
        self.assertTrue(all(row["singleton_matches_stored_truth"] for row in checked))
        self.assertTrue(all(row["nonzero_relevant_quotients"] == 11 for row in checked))


if __name__ == "__main__":
    unittest.main()
