# Archived experimental results

The archived release contains 1,050 runs in 35 configurations. Each configuration
uses three fixed secret values and ten inputs per secret; the 30 runs should not
be read as 30 independently sampled secrets. This report reaggregates the existing
records and compares saved outputs with saved truth. It does not rerun attacks,
lattice reduction, decoding, or the mathematical verification checks.

## Lowest tested leakage with 30/30 singleton recoveries

Here `b` is the number of unrevealed bits in a 384-bit field, and the revealed
fraction is `(384-b)/384`. Each row has 30 completed runs whose saved output is
exactly the singleton containing the saved true secret. Times are seconds; each
time cell is the median / maximum of the named original field.

| Problem / query model | b | Revealed | Full input → lattice input | Rank / ambient dimension | Parent time | Worker time |
|---|---:|---:|---:|---:|---:|---:|
| MIHNP / independent | 173 | 54.947917% | 30 → 12 | 78 / 78 | 6.357750 / 19.983134 | 5.713649 / 19.650027 |
| EC-HNP / paired + central | 124 | 67.708333% | 17 → 8 | 36 / 36 | 1.714786 / 2.793877 | 1.239952 / 2.297583 |
| EC-HNP / independent | 124 | 67.708333% | 17 → 8 | 36 / 36 | 1.711552 / 3.070154 | 1.172469 / 2.721744 |
| HNPHM / independent | 123 | 67.968750% | 16 pairs → all 16 pairs | 33 / 49 | 1.115435 / 9.193445 | 0.654538 / 4.178872 |

These are endpoints among the tested configurations under the recorded input and
resource settings. They do not establish a universal leakage threshold or an
asymptotic running-time guarantee. The archive contains no ours MIHNP run with
`b > 173`.

All 240 completed runs of our method, including the rows above, record
`nodes_visited=1` and `decoder_calls=1`. This is a root-node visit and one completed
decoder call. A final grid-leaf count was not saved, so this observation is not
reported as a first-grid-cell measurement. The additional readout and constraint
verification flags in the archive remain recorded evidence; this reaggregation
does not recertify them.

## Interrupted runs and baseline outcomes

The 150 `CENSORED_WORK` runs all record `nodes_visited=4096` and `secrets=null`:
30 runs each for EC-HNP paired at `b=125,126`, EC-HNP independent at `b=125,126`,
and HNPHM independent at `b=124`. They stopped at a work cap with no completed
output list. They are not mathematical impossibility results or recorded timeouts.
Their times and counters are summarized only in the `interrupted.*` columns.
The counter measures tree nodes, including internal regions, rather than final
grid cells.

| Archived prior implementation | Query model | Full lattice input | Rank / ambient dimension | Singleton-true counts, by b |
|---|---|---:|---:|---|
| MIHNP `lssw12` | independent | 30 observations | 60 / 60 | 122–124: 30/30 each; 125: 25/30; 126: 14/30; 127, 172, 173: 0/30 each |
| EC-HNP `shani17` | paired + central | 17 observations | 27 / 27 | 58–60: 30/30 each; 61–62 and 124: 0/30 each |
| HNPHM `hgns03` | independent | 16 pairs | 17 / 17 | 90: 30/30; 91: 29/30; 92: 20/30; 93: 6/30; 94: 2/30; 95: 1/30; 96 and 123: 0/30 each |

There are no independent-query EC-HNP prior runs in this archive. HNPHM at `b=90`
has a 76.5625% revealed fraction and is included above. These are results of the
stored implementations at the fixed input budget, rather than a claim about every
finite-size sufficient condition of the cited methods.

Across all runs the recorded statuses are `COMPLETE=240`, `CENSORED_WORK=150`,
`RETURNED=307`, and `NO_CANDIDATE=353`. There are 900 finished runs, 547
singleton-true outputs, and 353 other finished runs. A prior `NO_CANDIDATE` is an
algorithm termination, not a complete enumeration certificate. No `TIMEOUT` or
`ERROR` appears in these 1,050 records.

## Inputs, subset rule, and units

Our MIHNP and EC-HNP lattices use the fixed first 12 and first 8 observations,
respectively, and the released worker checks returned candidates against the full
input. HNPHM uses all 16 pairs, corresponding to 33 scalar reports including the
secret report. The prior implementations use the full inputs shown above; their
lattice input counts are not replaced with our subset sizes. For paired EC-HNP,
the first eight observations consist of the central observation, three complete
pairs, and one observation of the next pair. They are not eight independent
queries or four complete pairs.

The prefix rule is visible in the released code and agrees with the saved source
input counts. Per-run selected-index logs and the original rationale or tuning
history for choosing 12 and 8 were not retained. The prefix is reported as a
reconstruction from that code, not as a separately logged run parameter.

`source_rank` and `ambient_dimension` are distinct recorded quantities. HNPHM
input counts are pairs, while its scalar report count is `2N+1`. Times are seconds;
node and decoder counters are counts. `NA` in the CSVs means absent or
inapplicable and is never filled with zero. Each statistic has an `n_recorded`
column. `finished.*` includes all algorithm terminations, `singleton.*` includes
only correct singleton outputs, and `interrupted.*` includes only unfinished
work-capped runs.

## Timing and retained provenance

The two elapsed-time fields have different boundaries. In the released runner,
`parent_elapsed_seconds` starts before temporary request creation and subprocess
launch and ends after the worker result is read. It includes process startup,
imports, worker computation, and result I/O, but excludes the later parent-side
truth comparison and post-run diagnostics. In the released worker,
`worker_elapsed_seconds` starts after imports, request reading, and initial setup,
and ends before final result serialization. It includes preparation, recovery,
full-input candidate checking, and intermediate progress I/O. These times are
not added together, and their difference is not a measured phase breakdown.

The released runner has a 600-second default timeout. The per-run launch arguments
were not retained, so the default cannot be independently confirmed for every
record. The released worker uses a 4,096-node cap, consistent with all 150 observed
work-capped stops. Exact execution-time CPU/RAM, dependency and backend versions,
phase times, input-generation seeds, and a binding from each timed run to exact
code bytes were not retained. The released snapshot does not fully explain the
generation of every extra saved verification field. No current environment or
commit is substituted for the missing historical provenance.

## Recreate the reports without running experiments

From the repository root, use Python 3 and only its standard library:

```sh
python3 tools/summarize_archived_results.py --repo-root . --out results/reporting
```

The script reads the existing `results/trials.jsonl`, `config/cases.jsonl`, and
`config/truth.jsonl` without modifying them. It requires matching unique IDs,
checks the three-secret/ten-input cohorts and all saved output-size, truth-inclusion,
and singleton flags, and checks the archived status and counter totals. It writes
only the following derived files; no saved secrets or public inputs are copied
into them:

- [`configurations.csv`](../results/reporting/configurations.csv): all 35 configurations,
  with finished, correct-singleton, and interrupted statistics separated.
- [`endpoints.csv`](../results/reporting/endpoints.csv): the four rows in the first table.
- [`provenance.json`](../results/reporting/provenance.json): SHA-256 hashes of the three
  source files, validation counts, units, and scope.
