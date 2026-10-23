"""Aggregate experiment measurements by problem, method, and bit width."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
GROUP = ('family', 'N', 'observation_model', 'method', 'unknown_bits')


def read_lines(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'results/trials.jsonl')
    parser.add_argument('--cases', type=Path, default=ROOT / 'config/cases.jsonl')
    parser.add_argument('--out', type=Path, default=ROOT / 'results/summary.csv')
    args = parser.parse_args()
    cases = read_lines(args.cases)
    records = read_lines(args.input)
    by_id = {r['id']: r for r in records}
    if len(by_id) != len(records):
        raise ValueError('Repeated trial ID in results')
    if set(by_id) - {c['id'] for c in cases}:
        raise ValueError('Results contain trial IDs absent from the input cases')
    groups = defaultdict(list)
    for case in cases:
        groups[tuple(case[k] for k in GROUP)].append(case)
    rows = []
    for key, planned in sorted(groups.items()):
        measured = [by_id[c['id']] for c in planned if c['id'] in by_id]
        row = dict(zip(GROUP, key))
        field_bits = planned[0]['public']['p'].bit_length()
        row.update(known_bits=field_bits-row['unknown_bits'],
                   known_percentage=100*(field_bits-row['unknown_bits'])/field_bits,
                   planned=len(planned), executed=len(measured),
                   singleton_true=sum(r.get('singleton_true') is True for r in measured),
                   recovered_true_secret=sum(r.get('contains_true_secret') is True for r in measured),
                   complete=sum(r['status'] == 'COMPLETE' for r in measured),
                   output_failure=sum(r.get('secrets') is not None
                                      and r.get('contains_true_secret') is False
                                      for r in measured),
                   node_censored=sum(r['status'].startswith('CENSORED') for r in measured),
                   timeout=sum(r['status'] == 'TIMEOUT' for r in measured),
                   error=sum(r['status'] == 'ERROR' for r in measured))
        times = [r['parent_elapsed_seconds'] for r in measured if r.get('singleton_true')]
        row['successful_median_seconds'] = median(times) if times else ''
        for name, field in [('actual_source_rank', 'source_rank'),
                            ('source_observations', 'source_observations_used')]:
            values = [r[field] for r in measured if r.get(field) is not None]
            row[name+'_min'] = min(values) if values else ''
            row[name+'_max'] = max(values) if values else ''
        rows.append(row)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'{len(records)} runs, {len(rows)} configurations: {args.out}')


if __name__ == '__main__':
    main()
