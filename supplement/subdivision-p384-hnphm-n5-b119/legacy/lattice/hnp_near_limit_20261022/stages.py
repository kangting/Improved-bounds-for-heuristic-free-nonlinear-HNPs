"""Freeze paired control pilots, all270main inputs, or dimensional auxiliaries."""
import argparse
from copy import deepcopy
import json
from pathlib import Path

from common import ROOT, digest, read, save, sha, uniform_below, utc
from domain import P, M, PARAMETERS, validate
from plan import MAIN, AUXILIARY, exact_bits
from solver import core


def fixed_secret(family, material):
    if family == "echnp":
        ec = core.curve({"p": P, **PARAMETERS})
        scalar = 1 + uniform_below(M - 1, material)
        return core.encode_point(scalar * core.point(ec, PARAMETERS["generator"]))
    return uniform_below(P if family == "mihnp" else P - 1, material) + (family == "hnphm")


def make_record(family, n, bits, seed, secret, metadata):
    public, opening = core.generate(family, P, n, bits, seed, secret,
                                    PARAMETERS if family == "echnp" else None)
    meta = {"problem": family, "method": "flatter-ntl", "implementation": "ours",
            "N": n, "B": 384, "unknown_bits": bits, "X": 1 << bits,
            "known_fraction": (384 - bits) / 384,
            "actual_query_count": n if family != "hnphm" else 2 * n + 1,
            "expected_source_rank": 2 * n + 1 if family == "hnphm" else n * (n + 1) // 2,
            "expected_ambient_dimension": 3 * n + 1 if family == "hnphm" else n * (n + 1) // 2,
            **metadata}
    identifier = digest([public, meta])[:24]
    return ({"id": identifier, "public": public, "public_sha256": digest(public), "metadata": meta},
            {"id": identifier, "metadata": meta, **opening})


def controls(parent):
    publics = [json.loads(line) for line in (parent / "public.jsonl").read_text().splitlines()]
    truth = {row["id"]: row for row in map(json.loads, (parent / "private/truth.jsonl").read_text().splitlines())}
    records, openings = [], []
    for original in publics:
        old = truth[original["id"]]
        for offset in [-2, -4]:
            metadata = {**original["metadata"], "campaign": "ours_paired_control_pilot",
                        "offset": offset, "parent_pilot_id": original["id"],
                        "paired_observation_id": original["id"]}
            # Width-dependent metadata is set by make_record, not inherited.
            for key in ["unknown_bits", "X", "known_fraction"]:
                metadata.pop(key)
            record, opening = make_record(original["public"]["family"], original["public"]["N"],
                  original["metadata"]["unknown_bits"] + offset, old["observation_seed"], old["secret"], metadata)
            for key in ["shifts", "queries"]:
                if key in original["public"]:
                    assert record["public"][key] == original["public"][key]
            old_values = [None if e is None else r + e for r, e in zip(original["public"]["reports"], old["errors"])]
            new_values = [None if e is None else r + e for r, e in zip(record["public"]["reports"], opening["errors"])]
            assert old_values == new_values
            records.append(record)
            openings.append(opening)
    return records, openings, {"phase": "paired_control_pilot", "seeds_reused_only_with_original_pilots": True}


def fresh(stage):
    assert stage in {"main", "auxiliary"}
    master = f"hnp-p384-prompt-{stage}-20261022-v1"
    records, openings = [], []
    settings = MAIN if stage == "main" else AUXILIARY
    for family, n in settings:
        b0 = exact_bits(P, core.finite_exponent(family, n))
        for secret_index in range(3):
            secret = fixed_secret(family, [master, "fixed-secret", family, n, secret_index])
            for draw in range(10 if stage == "main" else 1):
                seed = int(digest([master, "observations", family, n, secret_index, draw]), 16)
                pair = digest([master, "paired-observations", family, n, secret_index, draw])[:24]
                values = []
                for offset in [0, -2, -4] if stage == "main" else [0]:
                    record, opening = make_record(family, n, b0 + offset, seed, secret,
                        {"campaign": "ours_main" if stage == "main" else "ours_auxiliary_pilot",
                         "fixed_secret_index": secret_index, "draw_index": draw, "offset": offset,
                         "secret_identifier": f"{stage}-{family}-N{n}-secret{secret_index}",
                         "paired_observation_id": pair})
                    values.append([None if e is None else r + e for r, e in zip(record["public"]["reports"], opening["errors"])])
                    records.append(record)
                    openings.append(opening)
                assert all(value == values[0] for value in values)
    return records, openings, {"master": master, "phase": stage, "never_passed_to_solver": True}


def freeze(parent_path, out_path, stage):
    parent, out = Path(parent_path).resolve(), Path(out_path).resolve()
    old_manifest = read(parent / "manifest.json")
    assert sha(parent / "public.jsonl") == old_manifest["public_sha256"]
    assert sha(parent / "private/truth.jsonl") == old_manifest["truth_sha256"]
    # Validate fixed public domain; proven-primality evidence remains bound below.
    validate(prove=False)
    records, openings, seeding = controls(parent) if stage == "controls" else fresh(stage)
    assert len(records) == {"controls": 18, "main": 270, "auxiliary": 18}[stage]
    assert len({row["id"] for row in records}) == len(records)
    out.mkdir(parents=True, exist_ok=False)
    (out / "private").mkdir()
    (out / "jobs").mkdir()
    save(out / "domain.json", read(parent / "domain.json"))
    for path, rows in [(out / "public.jsonl", records), (out / "private/truth.jsonl", openings)]:
        with path.open("x") as handle:
            for row in rows:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
    manifest = deepcopy(old_manifest)
    manifest.update(created_utc=utc(), planned_trials=len(records), stage=stage,
         parent_manifest_sha256=sha(parent / "manifest.json"),
         public_sha256=sha(out / "public.jsonl"), truth_sha256=sha(out / "private/truth.jsonl"),
         domain_sha256=sha(out / "domain.json"), source_sha256={**old_manifest["source_sha256"],
               "stages.py": sha(__file__), "runner_unbounded.py": sha(ROOT / "runner_unbounded.py"),
               "runner_unbounded_v2.py": sha(ROOT / "runner_unbounded_v2.py")},
         effective_per_job_timeout_seconds=None, individual_timeout_user_amendment="2026-10-22",
         inherited_timeout_seconds_is_historical_only=True,
         shared_campaign_clock_sha256=sha(parent.parent / "campaign_clock.json"),
         phase_input_plan="All frozen records remain planned; resource skips/not-run must be explicit and never count as successes")
    save(out / "manifest.json", manifest)
    save(out / "private/seeding.json", seeding)
    print(json.dumps({"stage": stage, "frozen_inputs": len(records), "out": str(out)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--stage", choices=["controls", "main", "auxiliary"], required=True)
    args = parser.parse_args()
    freeze(args.parent, args.out, args.stage)
