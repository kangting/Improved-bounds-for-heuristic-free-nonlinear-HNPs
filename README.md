# HNP experiments

Code and results for integer-quotient recovery of MIHNP, EC-HNP, and HNPHM,
with lattice-based baselines. The implementation uses Python, python-flint,
fpylll, and [flatter](https://github.com/keeganryan/flatter) for lattice
reduction, and Matplotlib for the figures.

## Files

- `src/`: polynomial relations, tracked lattice reduction, search, decoders,
  and baseline methods.
- `config/`: public experiment inputs and separate evaluation secrets.
- `results/`: 1,050 main-run measurements and their summary.
- `fig/`: leakage and recovery plots in PDF format.
- `run.py`, `worker.py`: experiment runner and isolated solver process.
- `summarize.py`, `plot.py`: result aggregation and plotting.

## Setup

The experiments were run with Python 3.12 on Linux.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install flatter using its upstream installation instructions. Put its executable
on `PATH`, or specify it directly:

```bash
export HNP_FLATTER_BIN=/path/to/flatter
```

## Recreate the figures

Plotting the stored results requires only Matplotlib.

```bash
python summarize.py
python plot.py
```

The plots are written to `fig/leakage.pdf` and `fig/recovery.pdf`.

## Run experiments

Run one stored MIHNP input:

```bash
python run.py --family mihnp --method ours --unknown-bits 173 \
  --limit 1 --out results/test.jsonl
```

Run all stored inputs and summarize the new measurements:

```bash
python run.py --workers 4 --out results/new.jsonl
python summarize.py --input results/new.jsonl --out results/new_summary.csv
python plot.py --data results/new_summary.csv --out fig/new
```

Use `--help` for filtering options. The runner does not overwrite existing
measurement files.
Each run has a 600-second default limit; our search also has a 4,096-node limit.
An unfinished search returns no candidate list.

The stored comparison uses 384-bit fields and 30 runs per configuration, with
public observation counts 30, 17, and 16 for MIHNP, EC-HNP, and HNPHM.
Our MIHNP and EC-HNP lattices use the first 12 and 8 observations; every
candidate is checked against the full input. HNPHM uses all 16 pairs.
The main EC-HNP comparison uses paired queries, with independent-query runs
included separately. Figure endpoints are the lowest tested leakage levels
with 30/30 unique recoveries.
