"""Isolated public-only cell-tree worker; explicit work caps, no job timer."""
import argparse
import json
from pathlib import Path
import resource
from time import monotonic
import traceback

from common import read, atomic, utc, digest
from solver_cells import solve_cells


def main(args):
    public, out = read(args.public), Path(args.out)
    if public["family"] == "echnp":
        from domain import validate
        validate(prove=False)
    resource.setrlimit(resource.RLIMIT_AS, (args.address_space, args.address_space))
    launched = monotonic()
    origin = read(Path(args.reuse_preparation) / "progress.json")["source_start_monotonic"] if args.reuse_preparation else launched
    state = {}
    with (out / "events.jsonl").open("x", buffering=1) as log:
        def event(**kw):
            state.update(kw)
            row = {**state, "utc": utc(), "monotonic": monotonic(),
                   "elapsed_seconds": monotonic() - origin,
                   "new_worker_elapsed_seconds": monotonic() - launched}
            log.write(json.dumps(row, sort_keys=True) + "\n")
            atomic(out / "progress.json", row)
        event(phase="solver_started", source_start_monotonic=origin,
              new_worker_start_monotonic=launched, algorithm_arm="adaptive_cell_tree")
        try:
            result = solve_cells(public, args.backend, event, out / "certificate.json.gz",
                         args.maximum_nodes, not args.no_original_pruning, args.reuse_preparation,
                         args.maximum_targets)
            result.update(public_sha256=digest(public), solver_seconds=monotonic() - origin,
                          new_worker_seconds=monotonic() - launched,
                          original_preparation_included_in_wall_time=bool(args.reuse_preparation))
            atomic(out / "solver_result.json", result)
            event(phase="complete_output" if result["status"] == "COMPLETE" else "resource_limited_partial",
                  terminal=True, complete=result["status"] == "COMPLETE")
        except BaseException as exc:
            event(phase="solver_error", terminal=True, error=repr(exc), traceback=traceback.format_exc())
            raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--backend", default="flatter-ntl", choices=["flatter-ntl", "ntl-exact"])
    parser.add_argument("--address-space", type=int, default=8 * 1024**3)
    parser.add_argument("--maximum-nodes", type=int, default=10000)
    parser.add_argument("--maximum-targets", type=int, default=10000)
    parser.add_argument("--no-original-pruning", action="store_true")
    parser.add_argument("--reuse-preparation")
    args = parser.parse_args()
    main(args)
