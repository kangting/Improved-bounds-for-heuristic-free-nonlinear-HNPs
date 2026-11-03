"""Evaluation-only stage, called after all public interval/replay checks finish."""
import json
from pathlib import Path

from verify_public import IDS, P, X, load, require, sha, verify_integrity


def evaluate(root, public_artifacts):
    # This is the first stage allowed to read the truth bytes.
    verify_integrity(root, evaluation_only=True)
    def keyed(path):
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        require(len(rows) == 3 and {row["id"] for row in rows} == set(IDS), "Evaluation selection is not exactly these three trials")
        return {row["id"]: row for row in rows}
    truths = keyed(root / "evaluation/truth.jsonl")
    historical = keyed(root / "evaluation/trials.jsonl")
    summary = load(root / "evaluation/evaluator_metadata.json")
    require(summary["evaluator_sha256"] == sha(root / "legacy/lattice/hnp_near_limit_20261022/evaluate.py"), "Historical evaluator source hash mismatch")
    output = []
    for identifier in IDS:
        truth, row = truths[identifier], historical[identifier]
        artifact = public_artifacts[identifier]
        job = root / "records/trials" / identifier
        public, result = load(job / "public.json"), load(job / "result.json")
        solver = result["solver_result"]
        require(solver["outputs"] == [truth["secret"]], "Stored output is not the stored true singleton")
        require(row["singleton_true"] is True and row["contains_true_secret"] is True and row["all_outputs_original_constraint_valid"] is True, "Historical evaluator flags disagree")
        require(row["result_sha256"] == sha(job / "result.json"), "Historical evaluator result hash mismatch")
        errors, multipliers = truth["errors"], truth["multipliers"]
        require(len(errors) == 11 and len(multipliers) == 5, "Invalid truth dimensions")
        values = [truth["secret"]] + multipliers + [truth["secret"] * t % P for t in multipliers]
        require(all(0 < value < P for value in values), "Invalid nonzero synthetic truth")
        require([value - report for value, report in zip(values, public["reports"])] == errors, "Truth errors do not match public reports")
        require(all(0 <= error < X for error in errors), "Truth witness is outside public intervals")
        polys = artifact["basis_evidence"]["full_H_reconstructed_now"]
        quotients = {}
        for j in artifact["basis_evidence"]["relevant_rows"]:
            polynomial = polys[j]
            value = polynomial[0] + sum(c * x for c, x in zip(polynomial[6:], errors))
            value += sum(c * errors[0] * errors[i + 1] for i, c in enumerate(polynomial[1:6]))
            q, remainder = divmod(value, P)
            require(remainder == 0, "Truth does not satisfy a reduced congruence")
            quotients[str(j)] = q
        alpha = artifact["basis_evidence"]["alpha"]
        require(sum(alpha[int(j)] * q for j, q in quotients.items()) == errors[0], "Post-public truth readout failed")
        cached = load(root / "evaluation/quotient_diagnostics" / (identifier + ".json"))
        require(cached["certificate_sha256"] == sha(job / "certificate.json.gz"), "Historical quotient diagnostic certificate binding mismatch")
        require(cached["relevant_quotients"] == quotients, "Historical quotient diagnostic differs")
        nonzero = sum(q != 0 for q in quotients.values())
        require(row["nonzero_relevant_quotient_count"] == nonzero and row["has_nonzero_relevant_quotient"] == (nonzero > 0), "Historical quotient flags differ")
        output.append({"id": identifier, "singleton_matches_stored_truth": True,
                       "public_replay_finished_before_truth_read": True, "nonzero_relevant_quotients": nonzero})
    return output
