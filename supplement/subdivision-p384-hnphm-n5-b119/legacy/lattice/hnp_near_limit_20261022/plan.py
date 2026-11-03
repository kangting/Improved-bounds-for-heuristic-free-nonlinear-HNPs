"""Freeze exact-window P384 pilots and the full requested experiment inventory."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from common import ROOT, digest, sha, save, utc, uniform_below
from domain import P, M, PARAMETERS, validate

MAIN = [("mihnp", 30), ("echnp", 18), ("hnphm", 16)]
AUXILIARY = [("mihnp", 12), ("mihnp", 20), ("echnp", 8),
             ("echnp", 12), ("hnphm", 5), ("hnphm", 11)]


def exact_bits(p, exponent):
    numerator, denominator = exponent.numerator, exponent.denominator
    bits = p.bit_length() * numerator // denominator
    while 2 ** (bits * denominator) > p ** numerator:
        bits -= 1
    assert 2 ** (bits * denominator) <= p ** numerator < 2 ** ((bits + 1) * denominator)
    return bits


def inventory():
    from solver import core
    rows = []
    for family, n in MAIN:
        b0 = exact_bits(P, core.finite_exponent(family, n))
        for offset in [0, -2, -4, 1, 2]:
            rows.append({"campaign": "ours_main" if offset <= 0 else "ours_optional_boundary",
                         "family": family, "N": n, "unknown_bits": b0 + offset,
                         "offset": offset, "fixed_secrets": 3, "draws_per_secret": 10,
                         "planned_trials": 30, "default_enabled": offset <= 0})
    for family, n in AUXILIARY:
        rows.append({"campaign": "ours_auxiliary", "family": family, "N": n,
                     "unknown_bits": exact_bits(P, core.finite_exponent(family, n)),
                     "planned_pilots": 3, "trigger": "main setting too expensive",
                     "default_enabled": False})
    for family, n, base_bits, count in [("mihnp", 30, 124, 30),
                                      ("echnp", 17, 60, 35), ("hnphm", 80, 72, 161)]:
        for offset in [0, -2, 2]:
            rows.append({"campaign": "prior_own_model", "family": family, "N": n,
                         "unknown_bits": base_bits + offset, "planned_trials": 30,
                         "actual_field_observation_count": count,
                         "default_enabled": False, "trigger": "source algorithm verified"})
    return rows


def freeze_pilot(out, backend="flatter-ntl"):
    from solver import core, FLATTER, FLATTER_LIBRARY, ADAPTER_DIR
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / "private").mkdir()
    (out / "jobs").mkdir()
    domain = validate(prove=True)
    save(out / "domain.json", domain)
    master = "hnp-p384-prompt-pilot-20261022-v1"
    records, openings = [], []
    for trial in range(3):
        for family, n in MAIN:
            bits = exact_bits(P, core.finite_exponent(family, n))
            material = [master, "fixed-secret", family, trial]
            if family == "echnp":
                ec = core.curve({"p": P, **PARAMETERS})
                scalar = 1 + uniform_below(M - 1, material)
                secret = core.encode_point(scalar * core.point(ec, PARAMETERS["generator"]))
            else:
                secret = uniform_below(P if family == "mihnp" else P - 1, material)
                secret += family == "hnphm"
            seed = int(digest([master, "observations", family, n, trial]), 16)
            public, opening = core.generate(family, P, n, bits, seed, secret,
                                            PARAMETERS if family == "echnp" else None)
            meta = {"problem": family, "campaign": "ours_pilot", "method": backend,
                    "implementation": "ours",
                    "N": n, "B": 384, "unknown_bits": bits, "X": 1 << bits,
                    "known_fraction": (384 - bits) / 384, "trial": trial,
                    "secret_identifier": f"pilot-{family}-{trial}",
                    "actual_query_count": n if family != "hnphm" else 2 * n + 1}
            identifier = digest([public, meta])[:24]
            records.append({"id": identifier, "public": public,
                            "public_sha256": digest(public), "metadata": meta})
            openings.append({"id": identifier, "metadata": meta, **opening})
    for path, data in [(out / "public.jsonl", records), (out / "private/truth.jsonl", openings)]:
        with path.open("x") as handle:
            for row in data:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
    sources = ["common.py", "domain.py", "plan.py", "solver.py", "worker.py", "runner.py"]
    save(out / "manifest.json", {"created_utc": utc(), "planned_trials": len(records),
         "stage": "ours_primary_pilot", "backend": backend,
         "public_sha256": sha(out / "public.jsonl"), "truth_sha256": sha(out / "private/truth.jsonl"),
         "domain_sha256": sha(out / "domain.json"),
         "prompt_sha256": sha(ROOT.parent / "HNP_Near_Limit_Experiment_Prompt_KO.md"),
         "paper_sha256": sha(ROOT.parent / "STOC2027_HNPs (19).pdf"),
         "source_sha256": {name: sha(ROOT / name) for name in sources},
         "dependency_sha256": {str(ADAPTER_DIR / "model.py"): sha(ADAPTER_DIR / "model.py"),
                               str(ROOT.parent / "hnp_threshold_experiments/hnp.py"): sha(ROOT.parent / "hnp_threshold_experiments/hnp.py")},
         "flatter_executable_sha256": sha(FLATTER), "flatter_library_sha256": sha(FLATTER_LIBRARY),
         "timeout_seconds": 600, "campaign_budget_seconds": 86400,
         "peak_rss_limit_bytes": 4 * 1024**3, "address_space_limit_bytes": 8 * 1024**3,
         "maximum_workers": 2, "numerical_threads_per_worker": 1,
         "main_preflight_policy": "Primary b0 pilot first; if costly, evaluate paired b0-2/b0-4 pilots before gating any width. Retain all30planned draws per main width; resource-preflight skips are not mathematical failures.",
         "recovery_private_input": False, "post_completion_diagnostics_only": True,
         "inventory": inventory(), "optional_field_sizes": [256, 521],
         "common_input_comparison": "not planned: observation models differ; any future matched-access arm requires separate commitment"})
    save(out / "private/seeding.json", {"master": master, "phase": "pilot", "never_passed_to_solver": True})
    print(json.dumps({"frozen_pilot_inputs": len(records), "out": str(out)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("--backend", choices=["flatter-ntl", "ntl-exact"], default="flatter-ntl")
    args = parser.parse_args()
    freeze_pilot(args.out, args.backend)
