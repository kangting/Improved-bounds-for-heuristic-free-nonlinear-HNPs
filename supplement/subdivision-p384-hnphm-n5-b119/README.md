# Existing P-384 HNPHM subdivision example: N=5, b=119

This is a separate auxiliary supplement, not part of the 1,050-run main
benchmark or its N=16 HNPHM success/time statistics. It publishes all three
previously recorded inputs at this auxiliary setting, their recovered
original integer certificates, and an independently inspectable check-only
implementation. No new inputs, LLL reduction, attack experiment or runtime
measurement was performed to produce this supplement. The paper source is
not changed.

## Check without Sage or lattice software

From the repository root, using Python3.9+ and only its standard library:

```sh
python3 -B supplement/subdivision-p384-hnphm-n5-b119/check.py
python3 -B supplement/subdivision-p384-hnphm-n5-b119/check.py --public-only
python3 -B -m unittest discover -s supplement/subdivision-p384-hnphm-n5-b119/tests -p 'test_*.py'
python3 -B audit_results.py --check
```

The first command checks public material first, then runs the separate
evaluation-only stage. `--public-only` does not open anything in `evaluation/`,
including truth, seeds, errors and cached private diagnostics. Neither mode
imports the preserved legacy modules or calls an LLL routine or legacy
decoder. The checker verifies **stored** candidates; it does not run fresh
reconstruction branches. Tests include a guard rejecting any public-stage
attempt to open evaluation files, as well as certificate/trace tampering.

Expected per-trial output (long hashes and additional PASS lines omitted):

```text
7129544b15d539fb71e5e983: root_targets=1 verification_nodes=1 stored_decoder_calls=1 empty_target_prunes=0
210cec052c05c4c5e7ac8845: root_targets=175656115840599632221463450780895704 verification_nodes=5 stored_decoder_calls=1 empty_target_prunes=2
4b2bf95e230d9cfa2a6c1929: root_targets=1 verification_nodes=1 stored_decoder_calls=1 empty_target_prunes=0
PASS: evaluation-only stage 3/3 stored singletons, original witnesses and quotient diagnostics
```

`stored_decoder_calls` describes the original recorded call counter, checked
against the unique singleton-target requests in the verification replay.
It is not a claim that the new checker called a decoder. Failures produce
nonzero exit status. Default checks do not write files. Optionally:

```sh
python3 -B supplement/subdivision-p384-hnphm-n5-b119/check.py --public-only --write-replay work/subdivision-verification
```

The destination must be new. Those files are another copy of the **new
verification replay**, without timings, not newly recovered historical logs.

## Historical setting and results

The prime is `2^384 - 2^128 - 2^96 + 2^32 - 1`. There are five hidden
multiplier/product pairs and one secret approximation: eleven scalar
reports, rank11 and ambient dimension16. `X=2^119` is retained as the
historical radius/scale parameter, with the declared signed box `[-X,X]^11`
intersected at each node with the original public leading-bit intervals.
No planted error is used for those intersections.

The original backend is flatter preprocessing followed by exact tracked NTL
LLL at delta3/4, not the main benchmark's later FLINT backend. The historical
adaptive policy has a **10,000-node** limit, no individual wall timeout,
a separate bounded-target safety limit, and a shared 24-hour campaign clock.
These three runs finish far below the node limit with no target fallback.
The original campaign also used sampled process-group RSS and address-space
limits; the original settings remain in `records/campaign/`.

| Trial ID | Historical nodes | Historical target/decoder calls | Historical outcome |
|---|---:|---:|---|
|7129544b15d539fb71e5e983|1|1|COMPLETE, stored singleton true|
|210cec052c05c4c5e7ac8845|5|1|COMPLETE, stored singleton true|
|4b2bf95e230d9cfa2a6c1929|1|1|COMPLETE, stored singleton true|

For the middle case, the exact public root target interval is

```text
[64107237951601111634317879330044870,
 239763353792200743855781330110940573]
```

It contains `175656115840599632221463450780895704` integers (about
1.76e35). This is the intersection of the sign-aware interval image
`sum(alpha[j] * J[j])` with the first public coordinate interval. It is not
obtained by substituting the true errors. The relevant quotient intervals
include one nonsingleton interval, but ambiguity of individual quotients
alone is not used as proof of ambiguity of the weighted target.

## What was recovered unchanged

`records/trials/<id>/` contains each original:

- `public.json`, `certificate.json.gz`;
- `result.json`, `solver_result.json`, `progress.json`;
- `events.jsonl`, and a documented path-normalized copy of `process.json`.

All three certificates were recovered from the original working tree, not
recreated or newly reduced. Their gzip bytes match the hashes already saved
in the earlier `provenance/prior_audit/p384_auxiliary.csv` and the historical
post-completion quotient diagnostics. The original certificates store **U
and alpha**, not a separate reduced-basis R array. The public source B0 and
R=U*B0 are reconstructed now by exact integer multiplication. The new basis
arrays are clearly labelled `*_reconstructed_now`; they are not fabricated
original basis files.

`records/public.selected.jsonl` preserves the three original public-row
lines, including metadata and public digests. `evaluation/truth.jsonl` and
`evaluation/trials.jsonl` preserve only the three requested original lines.
The truth contains synthetic secret, errors, multipliers and seed, and is
evaluation-only. The three original cached quotient diagnostic files and
the historical evaluator metadata are also retained under `evaluation/`.
No other campaign's private openings are included.

## Verification stages and mathematical scope

1. Check publication SHA-256 values and the fifteen original campaign
   source/dependency bindings. The evaluator has a separate historical
   evaluator-summary binding, not a campaign-source binding.
2. Rebuild all original HN polynomials, including constants and target-last
   ordering, and their 11x16 scaled source lattice B0. Check a nonzero rank
   minor, integral U, `abs(det(U))=1`, full R=U*B0, alpha*U=e_target and the
   complete polynomial target identity. Rebuild relevant-row W and the
   recorded subdivision parameter.
3. Using only public inputs and that certificate, compute every relevant
   polynomial/quotient interval at the root and at each replay node. Apply
   signs of alpha before intersecting with the public first-coordinate
   interval. Reproduce the original balanced-grid boundaries, deterministic
   coefficient-sensitivity split rule, DFS order and pruning decisions.
4. Check the new replay's coverage and aggregate counters against saved
   records, and verify stored outputs against the replay's singleton-target
   readout and all original input intervals. Hidden-multiplier existence
   for a **stored candidate** is checked with exact logarithmic floor sums,
   not an interval scan or a new decoder/recovery call.
5. Only after all public checks complete, open the selected truth and
   historical evaluator data. Compare the saved singleton output, original
   witness and post-completion quotient/readout diagnostics. Truth cannot
   change the already computed intervals, splits or pruning decisions.

This does not run fresh LLL or independently certify its large rational GSO
quality. It checks exact algebraic identities of the saved certificate and
the interval traversal. The historical NTL LLL contract remains a library
contract, not a newly measured/independently proven reduction-quality claim.
No new timing or new sampling-distribution evidence is provided.

## Historical events versus the NEW verification replay

The original `events.jsonl` files contain phases and aggregate checkpoints.
They do **not** record each node's box, polynomial bounds or quotient bounds.
No missing detailed historical trace is claimed to have been recovered.

`verification_replay/<id>.json` is generated now from the saved public input,
certificate and preserved historical split rule. It includes B0/R/H,
per-node grid ranges and boxes, every processed relevant polynomial and
quotient interval, pre/post-intersection target bounds, split scores,
decisions, logical cell coverage and pending frontier. It contains no
elapsed-time measurement. The checker recomputes these files and compares
them byte-for-byte, in addition to checking their hashes.

For the five-node input its decision tree is:

```text
root: SPLIT coordinate0
  L: EMPTY_TARGET_INTERVAL
  R: SPLIT coordinate0
    L: SINGLETON_TARGET
    R: EMPTY_TARGET_INTERVAL
```

This newly checked rule-consistent trace has five nodes, two empty-target
prunes, one unique singleton-target request and frontier zero, matching the
historical aggregates. All `17^11=34271896307633` logical grid cells are
covered by the subtree certificates. That number is not a visited-node count
or a runtime. This is a nontrivial subdivision **mechanism example** on a
fixed saved basis. It is not a general speedup measurement, a universal
impossibility result for other root-only methods, or evidence attributing
the main N=16 benchmark's improvements causally to branching.

## Provenance, normalized metadata, and legacy source

`manifest.json` lists every published payload file, its SHA-256, role and
available workspace-relative source/selection information. It excludes its
own self-hash; the Git commit identifies this manifest. These new publication
hashes identify current files and are not a retrospective signature of the
original machine state or every historical runtime event.

`provenance/recovery.json` records original and published hashes and exact
metadata normalization fields. Public inputs, certificates, numerical
results, progress and event bytes are unchanged. Local deployment paths in
the three process records and campaign metadata are replaced with explicit
placeholders; a machine boot ID is similarly normalized. No numerical
measurement, transform, quotient, target, counter or scientific setting is
changed by this metadata-only normalization.

The original campaign manifest's public/truth hashes cover **eighteen**
campaign records, not the selected three. Those full-stream commitments are
preserved as historical metadata but are not claimed validated by the
reduced publication streams. Selected lines have their own publication
hashes and source-line provenance. The historical evaluator summary also
describes the larger checkpoint, not an independent three-case campaign.

`legacy/lattice/` preserves all thirteen manifest-listed project sources and
the two sibling dependencies byte-for-byte, plus separately bound evaluator
support. The original paths/import layout is retained:

```text
worker_cells -> solver_cells -> solver -> hnp_table1_v19/model
                                      -> hnp_threshold_experiments/hnp
campaign_cells -> runner_unbounded -> runner
               -> runner_unbounded_v2
freeze_cells/stages -> plan/domain/solver
```

The legacy sources retain historical deployment literals (including the
flatter binary/library paths) to preserve their original hashes. They need
SageMath10.8/NTL/flatter for a historical execution and are **not** imported
by this checker. Do not execute the archival campaign runners to check this
supplement: they are not a turnkey relocated campaign, contain historical
supervision/termination code, and require external runtime/configuration
assets not needed for these checks. Optional MI-only binary and non-ours
prior branches were not used for these three HN jobs and are not supplied.
No binary, manuscript PDF or executable runtime image is bundled.

The earlier `README_KO.md`, `audit_existing.py` and `p384_auxiliary.csv` were
located and copied unchanged into `provenance/prior_audit/`. They are
archival evidence of the previous wider local audit, not this supplement's
standalone entry point; the old script also refers to out-of-scope legacy
61-bit data. Use `check.py` for the self-contained public verification.

All source hash matches, mathematical rechecks and replay consistency
claims have the scopes above. Missing original per-node logging, original
runtime behavior beyond the retained records, and independent large-basis
LLL-quality certification remain unverified. The main benchmark files and
their success/time totals are not modified or augmented by these cases.
