"""Recover from one public input in an isolated process."""

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
sys.path.insert(0, str(ROOT / "src"))

from baselines_native import prior
from hnp_native import prepare, recover_pruned
from subinstance import public_subset, verify_candidates


def write_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, separators=(",", ":")))
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    request = json.loads(args.request.read_text())
    public, method = request["public"], request["method"]
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    metadata = {}

    def phase(name, data):
        write_json(out / "progress.json", {**metadata, "last_phase": name})

    try:
        if method == "prior":
            phase("prior", {})
            result = prior(public)
            metadata = {"source_rank": result.get("dimension"),
                        "source_observations_used": public["N"]}
        elif method == "ours":
            if public["family"] == "echnp" and None in public["reports"]:
                result = recover_pruned(public)
                metadata = {"source_rank": None,
                            "source_observations_used": public["N"]}
            else:
                active = (public_subset(public)
                          if public["family"] in ("mihnp", "echnp") else public)
                prepared = prepare(active, phase=phase)
                metadata = {
                    "source_rank": prepared["basis"].nrows(),
                    "source_observations_used": active["N"],
                }
                phase("search", {})
                result = recover_pruned(
                    active, prepared, phase=phase,
                    progress=lambda data: phase("search", data), max_nodes=4096)
                if result["status"] == "COMPLETE":
                    result["secrets"] = verify_candidates(public, result["secrets"])
                    # Used only by the evaluation parent after recovery; this
                    # temporary data is omitted from the saved trial results.
                    count = active["N"]
                    indices = (list(range(2 * count + 1))
                               if public["family"] == "hnphm" else list(range(count)))
                    result["_postrun"] = {
                        "error_indices": indices,
                        "edges": prepared["sources"].edges,
                        "alpha": prepared["alpha"],
                        "relevant": [
                            {"index": j,
                             "constant": prepared["polynomials"][j].constant,
                             "linear": prepared["polynomials"][j].linear,
                             "quadratic": prepared["polynomials"][j].quadratic}
                            for j in prepared["relevant"]
                        ],
                    }
                else:
                    result["secrets"] = None
        else:
            raise ValueError(f"Unknown method: {method}")

        # Keep the run output small: solver detail remains in the core API.
        kept = {key: result[key] for key in
                ("status", "secrets", "nodes_visited", "decoder_calls", "_postrun")
                if key in result}
        write_json(out / "result.json", {
            **metadata, **kept,
            "worker_elapsed_seconds": time.perf_counter() - started,
        })
    except Exception as exc:
        write_json(out / "result.json", {
            **metadata, "status": "ERROR", "secrets": None,
            "error": f"{type(exc).__name__}: {exc}",
            "worker_elapsed_seconds": time.perf_counter() - started,
        })
        raise


if __name__ == "__main__":
    main()
