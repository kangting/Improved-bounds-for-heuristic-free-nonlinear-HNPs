# Stored experiment records and post-hoc aggregation

The supplement in `results/audit/` describes only the 1,050 measurements
already stored in this repository. No earlier campaign, alternative modulus,
new attack, or saved-basis replay is added to these denominators. Input,
truth, measurement, original summary and figure files remain unchanged.

## Files and checks

| File | Meaning |
|---|---|
| `runs.csv` | Projection of original scalar fields, plus explicitly prefixed `audit.*` interpretations; saved secret arrays are not duplicated |
| `summary.csv` | Configuration-level completion, outcomes, work and separate timing statistics |
| `cohorts.csv` | Secret-index counts, trial indices and fixed-secret checks |
| `field_coverage.csv` | Presence/non-null/missing counts for original fields |
| `provenance.json` | Relative-path input hashes, current release-source hashes and aggregation-script hash |

`python3 audit_results.py --check` regenerates these five files in memory
and compares their bytes. It imports only the Python standard library.
Stored outputs are compared with stored evaluation secrets by equality;
this is not a new recovery, an independent constraint check, a quotient
certificate audit, or an exhaustive proof of large-field list completeness.
The already public `config/truth.jsonl` remains evaluation-only data, not a
solver input. No additional private openings are included.

Original scalar field names and meanings are preserved. Additional columns
start with `audit.`. `NA` denotes missing or inapplicable information; zero,
false, and an empty list remain distinct. The original JSONL is the authority
for fields not projected into the CSV.

## Configuration and selection

All 35 stored configurations contain 30 runs: three distinct fixed evaluation
secrets, each with trial indices 0 through 9. These checks do not establish
sampling independence from a missing PRNG seed or turn paired comparisons
into additional independent draws.

| Family | Full input | Our source | Our rank / ambient | Prior source and rank / ambient |
|---|---|---|---|---|
| MIHNP | 30 observations | first 12 | 78 / 78 | all 30; 60 / 60 |
| EC-HNP | 17 queries | first 8 | 36 / 36 | all 17 paired; 27 / 27 |
| HNPHM | 16 multiplier/product pairs and secret approximation | all pairs | 33 / 49 | all pairs; 17 / 17 |

The source counts are stored per run. The prefix indices in the new CSV are
reconstructed from the release's `src/subinstance.py` and `worker.py`, not
from a per-run indices log. Candidates from a completed subset search are
checked against the full public input. An incomplete subset search is not
returned as a complete list. The HN5 value in the subset helper is not used
by this worker: it calls that helper only for MIHNP and EC-HNP.

For paired EC queries the release generator uses the order central `3R0`,
then `+Q1,-Q1,+Q2,-Q2,...` around that central query. Its first eight queries
contain three complete pairs and one member of the next pair, besides the
central query. They are not four complete pairs or eight iid queries.
Independent-query EC runs are tabulated separately. There is no independent
EC prior run in this dataset; this absence is not a zero success rate.

No retained rationale explaining why 12 and 8 were selected was found in
these release records. The constants and their use are documented, but a
performance tradeoff or an optimal-subset claim must not be inferred.
The published dataset has fixed source sizes and a bit-width sweep. It does
not contain a same-input subset-size/baseline-size tuning history. This does
not establish that no unrecorded exploration ever occurred. The existence
of the current case list alone is not evidence of preregistration.

The original `plot.py` chooses, among recorded widths, the largest `b` with
30/30 singleton recoveries and no unfinished run. This is the lowest tested
all-success leakage for the fixed setup, not a universal threshold or the
best result over every possible source size.

## Outcomes and resource limits

The records contain 240 `COMPLETE`, 150 `CENSORED_WORK`, 307 `RETURNED`, and
353 `NO_CANDIDATE` outcomes. Thus 900 algorithms finished, 547 returned
exactly the stored true secret, and 353 finished without a candidate.
`RETURNED` and `NO_CANDIDATE` are prior single-answer algorithm outcomes,
not our complete-list contract. No `TIMEOUT` or `ERROR` appears in these
1,050 records. Empty returned lists may pass a vacuous output-validity flag;
that flag is not a recovery success metric.

The four requested all-success points are:

| Family / query model | Unrevealed bits `b` | Revealed (%) | Singleton | Parent median / max (s) |
|---|---:|---:|---:|---:|
| MIHNP / independent | 173 | 54.947917 | 30/30 | 6.357750 / 19.983134 |
| EC-HNP / paired | 124 | 67.708333 | 30/30 | 1.714786 / 2.793877 |
| EC-HNP / independent | 124 | 67.708333 | 30/30 | 1.711552 / 3.070154 |
| HNPHM / independent | 123 | 67.968750 | 30/30 | 1.115435 / 9.193445 |

EC `b=125,126` in both models and HN `b=124` account for the 150 censored
runs. Each records 4,096 visited nodes. This is work censoring, not a proof
of impossible recovery or a timeout. No partial-candidate arrays were
retained. Our MI `b=172` also has 30/30 recovery; no larger-than-173 our-MI
width is present here.

The prior functions are `lssw12`, `shani17`, and `hgns03`. Their results at
all tested widths remain in the supplementary summary, including HN prior
`b=90` with 30/30 at 76.5625% revealed, and partial prior success at wider
widths. A matched-budget instantiation is not automatically an instance of
all sufficient finite-size conditions in the prior paper. In particular,
HN16 is not the HGNS03 80-pair prescription for epsilon=.01.

The release runner has a 600-second default timeout for both methods. There
is no per-run launch manifest proving whether that default was overridden.
Our worker specifies `max_nodes=4096`; the measured censoring counts match
that constant. These are tree nodes, including internal/rejected regions,
not necessarily final grid cells. Priors do not use our node cap.

## Counter and timing semantics

All 240 completed our-method records have `nodes_visited=1` and
`decoder_calls=1`. These are root-region completions, not independently
recorded first-grid-cell completions. `final_grid_leaves_visited` is not
retained by the worker's output projection.

The current `recover_pruned.propose` deduplicates targets before fully
decoding them. Under that code, `decoder_calls` equals distinct fully tested
targets. A separate tested-target count or target set is not retained,
and exact timed-code provenance is not established, so the independently
recorded distinct-target and first-grid-cell columns remain `NA`. A separate
conditional code-based interpretation is clearly labelled in `runs.csv`.

- `parent_elapsed_seconds` starts before the temporary directory/request is
  created and ends after the child exits and its progress/result files are
  read and merged. It includes child startup/imports, request/output I/O,
  computation and teardown. It excludes subsequent parent truth comparison,
  quotient diagnostics, temporary-directory cleanup and final result-stream
  writing. Campaign queue time is not a per-run parent measurement.
- `worker_elapsed_seconds` starts after imports, argument/request parsing
  and output-directory creation. It includes source construction, reduction,
  exact transformation/readout checks, search, full-input candidate filtering,
  intermediate progress I/O and preparation of post-run coefficient data.
  It is sampled before final result serialization/writing and excludes parent
  evaluation and post-run quotient diagnostics.

These ranges describe the released code; the elapsed-time difference must
not be labelled pure startup overhead. The source API produces component
times and backend details, but the worker does not retain them and its phase
callback discards the data argument. Source/reduction/search/verification
times therefore cannot be recovered by subtracting totals. They remain `NA`.
New instrumentation would be required; no such experiment was run here.

The summary separates `finished_*`, `singleton_*` and `interrupted_*` times
and counters, each with an explicit recorded-value count, median and max.
Finished-only timings are conditional, not a full-input completion-time
average when a configuration contains censoring. Missing prior counters
are not replaced by zero.

## Backend and provenance limits

Our source and the priors use the shared `lattice_native.tracked_lll` route:
uniform gcd normalization, flatter `-delta 0.99 -v`, exact recovery/checking
of the full-ambient row transformation, then a finishing LLL. The finisher
tries FLINT exact-Gram delta=.99/eta=.5, approximate-Gram .99/.5, exact-Gram
.75/.5, and finally the Fraction exact fallback .75/.5. The first accepted
valid transformation ends this ordered retry path; this is not best-of
subset tuning. The chosen path is not saved per trial. Prior SVP uses an
LLL first-vector approximation; HGNS uses LLL/Babai, not exact CVP.

README/requirements specify Python3.12/Linux and pin python-flint0.9.0,
fpylll0.6.4 and cysignals1.12.6, but they are not a historical installation
receipt. The worker sets OMP/OpenBLAS threads to1; the flatter subprocess
also sets MKL/BLIS threads to1. Actual historical CPU, RAM, kernel, package
patch versions, flatter hash/version, launch command and concurrent worker
count are not retained. The runner's workers=1 default and README's
workers=4 example do not establish the measured concurrency.

Seeds and timestamps are not present in the saved case/truth/trial records.
Secret/trial IDs are not PRNG seeds. Current machine details, local paths
and unrelated older campaigns are deliberately not substituted for missing
historical information or included in this public supplement.

Stored results include `ambient_dimension`, `true_relevant_quotients_all_zero`,
`postrun_true_readout_verified`, and `outputs_original_constraints_valid`,
but the current worker/runner do not explain all their creation/preservation
paths. The records and their flags are preserved, not silently rewritten.
All completed our-method records retain nonzero quotient/readout indicators;
those diagnostics were not recomputed by this JSON-only audit. The release
source hash inventory describes available code, not a certificate that those
exact bytes produced every timed record. The publication base commit in
`provenance.json` has the same limited meaning.
