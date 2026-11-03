"""Check this auxiliary supplement, without running legacy attacks or LLL."""
import argparse
from pathlib import Path

from verify_public import IDS, VerificationError, check_public, json_bytes

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--public-only", action="store_true", help="Do not read any evaluation/truth file")
    parser.add_argument("--write-replay", type=Path, help="Optional NEW directory for a fresh verification copy; no timings")
    args = parser.parse_args()
    try:
        root = args.root.resolve()
        reports = check_public(root)
        print("PASS: public file hashes; 15 campaign source/dependency bindings; exact B0/U/R/H and target identities")
        for identifier in IDS:
            replay = reports[identifier]["verification_replay"]
            print(f"{identifier}: root_targets={replay['nodes'][0]['final_target_cardinality']} "
                  f"verification_nodes={replay['visited_nodes']} "
                  f"stored_decoder_calls={len(replay['distinct_singleton_targets'])} "
                  f"empty_target_prunes={replay['empty_target_interval_subtrees']}")
        print("PASS: new verification replay matches stored aggregate counters; no historical detailed trace is claimed")
        if args.write_replay is not None:
            args.write_replay.mkdir(parents=True, exist_ok=False)
            for identifier, artifact in reports.items():
                with (args.write_replay / (identifier + ".json")).open("xb") as handle:
                    handle.write(json_bytes(artifact))
            print("Wrote a NEW verification-replay copy, without runtime measurements")
        if args.public_only:
            print("PASS: public-only mode; evaluation/truth was not read")
        else:
            from evaluate_saved import evaluate
            checked = evaluate(root, reports)
            print(f"PASS: evaluation-only stage {len(checked)}/3 stored singletons, original witnesses and quotient diagnostics")
        print("No new input, LLL, legacy decoder, attack timing, or legacy runner was executed.")
    except (VerificationError, OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(1, f"FAIL: {type(exc).__name__}: {exc}\n")


if __name__ == "__main__":
    main()
