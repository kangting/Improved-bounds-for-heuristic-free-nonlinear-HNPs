"""One public-only solver, isolated by its parent's resource-bounded group."""
import argparse
import json
from pathlib import Path
import resource
from time import monotonic
import traceback

from common import read, atomic, utc, digest
from solver import solve


def work(public_path, out, backend, address_space, implementation="ours"):
    out = Path(out)
    public = read(public_path)
    if public["family"] == "echnp":
        from domain import validate
        validate(prove=False)  # Fixed public domain already proven at commitment.
    resource.setrlimit(resource.RLIMIT_AS, (address_space, address_space))
    started = monotonic()
    with (out / "events.jsonl").open("x", buffering=1) as log:
        state = {}
        def event(**kw):
            state.update(kw)
            row = {**state, "utc": utc(), "monotonic": monotonic(), "elapsed_seconds": monotonic() - started}
            log.write(json.dumps(row, sort_keys=True) + "\n")
            atomic(out / "progress.json", row)
        event(phase="solver_started", source_start_monotonic=started)
        try:
            if implementation == "ours":
                result = solve(public, backend, event, out / "certificate.json.gz")
            else:
                from priors import solve_prior
                result = solve_prior(implementation, public, backend, event)
            result.update(public_sha256=digest(public), solver_seconds=monotonic() - started)
            atomic(out / "solver_result.json", result)
            event(phase="complete_output", complete=True)
        except BaseException as exc:
            event(phase="solver_error", error=repr(exc), traceback=traceback.format_exc())
            raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--backend", choices=["flatter-ntl", "ntl-exact"], required=True)
    parser.add_argument("--address-space", type=int, default=8 * 1024**3)
    parser.add_argument("--implementation", choices=["ours", "lssw12", "shani17", "hgns03"], default="ours")
    args = parser.parse_args()
    work(args.public, args.out, args.backend, args.address_space, args.implementation)
