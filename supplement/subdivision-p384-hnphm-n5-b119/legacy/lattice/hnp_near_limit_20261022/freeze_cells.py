"""Commit an explicitly paired adaptive-cell arm without altering reference runs."""
import argparse
from copy import deepcopy
import json
from pathlib import Path

from common import ROOT, read, save, sha, digest, utc


def freeze(source_path, out_path):
    source, out = Path(source_path).resolve(), Path(out_path).resolve()
    old = read(source / "manifest.json")
    assert sha(source / "public.jsonl") == old["public_sha256"]
    assert sha(source / "private/truth.jsonl") == old["truth_sha256"]
    records = [json.loads(x) for x in (source / "public.jsonl").read_text().splitlines()]
    truth = {r["id"]: r for r in map(json.loads, (source / "private/truth.jsonl").read_text().splitlines())}
    out.mkdir(parents=True, exist_ok=False)
    (out / "jobs").mkdir()
    (out / "private").mkdir()
    new_public, new_truth = [], []
    for record in records:
        metadata = {**record["metadata"], "algorithm_arm": "adaptive_cell_tree",
                    "reference_record_id": record["id"], "reference_run": source.name,
                    "literal_algorithm1": False}
        identifier = digest([record["public"], metadata])[:24]
        new_public.append({**record, "id": identifier, "metadata": metadata})
        new_truth.append({**truth[record["id"]], "id": identifier, "metadata": metadata})
    for path, rows in [(out / "public.jsonl", new_public), (out / "private/truth.jsonl", new_truth)]:
        with path.open("x") as handle:
            for row in rows:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
    save(out / "domain.json", read(source / "domain.json"))
    save(out / "private/seeding.json", {"paired_reference_run": str(source),
         "source_seeding": read(source / "private/seeding.json"), "never_passed_to_solver": True})
    manifest = deepcopy(old)
    files = ["solver_cells.py", "worker_cells.py", "campaign_cells.py", "freeze_cells.py",
             "runner_unbounded_v2.py", "runner_unbounded.py"]
    manifest.update(created_utc=utc(), stage="adaptive_" + old["stage"],
          algorithm_arm="adaptive_cell_tree", literal_algorithm1=False,
          parent_manifest_sha256=sha(source / "manifest.json"), paired_reference_run=str(source),
          public_sha256=sha(out / "public.jsonl"), truth_sha256=sha(out / "private/truth.jsonl"),
          domain_sha256=sha(out / "domain.json"),
          source_sha256={**old["source_sha256"], **{name: sha(ROOT / name) for name in files}},
          maximum_nodes=10000, maximum_target_fallback=10000,
          original_interval_pruning=True, effective_per_job_timeout_seconds=None,
          maximum_new_preparation_cost_preflight_seconds=3600,
          per_run_maximum_workers=2, global_numerical_slots_including_preserved_preparation=3,
          shared_cpu_leases=[16, 18], preserved_preparation_cpu=14,
          sampled_global_rss_bound_bytes=3 * 4 * 1024**3,
          user_amendment="Try C with pruning first; only small B1 fallback; skip excessively costly work",
          work_cap_is_not_completion=True,
          routing_policy="No up-front min(B1,C): sign-aware target-range pruning in adaptive C first; switch to untested targets ONLY after visiting exactlyT=|B1|nodes with frontier remaining and within the target resource cap. External node cap<T never triggers fallback. Main widths run30times only when their3matching pilots complete; otherwise retain30explicit preflight skips.")
    save(out / "manifest.json", manifest)
    print(json.dumps({"frozen_inputs": len(new_public), "out": str(out), "arm": "adaptive_cell_tree"}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    freeze(args.source, args.out)
