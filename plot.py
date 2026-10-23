"""Reproduce the figures from the measured matched-N summary, without solvers."""
from pathlib import Path
import argparse
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
PANELS = [('mihnp', 'independent', 30, 'MIHNP'),
          ('echnp', 'paired', 17, 'EC-HNP'),
          ('hnphm', 'independent', 16, 'HNPHM')]
COLORS = {'prior': '#D55E00', 'ours': '#0072B2'}
MARKERS = {'prior': 's', 'ours': 'o'}


def load_rows(path):
    rows = []
    for raw in csv.DictReader(path.open()):
        if raw.get('stage', 'main') != 'main':
            continue
        row = dict(raw)
        for key in ['N', 'unknown_bits', 'planned', 'executed', 'singleton_true',
                    'node_censored', 'timeout', 'error']:
            row[key] = int(row[key])
        row['known_percentage'] = float(row['known_percentage'])
        assert 0 <= row['singleton_true'] <= row['executed'] <= row['planned']
        assert abs(row['known_percentage'] - 100 * (384-row['unknown_bits']) / 384) < 1e-10
        rows.append(row)
    return rows


def group(rows, family, model, method):
    selected = [r for r in rows if r['family'] == family
                and r['observation_model'] == model and r['method'] == method]
    assert selected
    return sorted(selected, key=lambda r: r['unknown_bits'])


def save(fig, output, name):
    for suffix in ['pdf']:
        fig.savefig(output / f'{name}.{suffix}', dpi=360, bbox_inches='tight',
                    facecolor='white', metadata={'Creator': 'Matplotlib'})
    plt.close(fig)


def legend_handles(include_budget=False):
    result = [Line2D([0], [0], marker=MARKERS[m], color='none',
                     markerfacecolor=COLORS[m], markeredgecolor='white',
                     markersize=7.5, label=label)
              for m, label in [('prior', 'Prior'), ('ours', 'Ours')]]
    if include_budget:
        result.append(Line2D([0], [0], marker='x', color='none', markeredgecolor='#777777',
                             markersize=6, label='Unfinished'))
    return result


def leakage_plot(rows, output):
    fig, ax = plt.subplots(figsize=(7.2, 2.85))
    fig.subplots_adjust(left=.25, right=.97, bottom=.24, top=.82)
    exported = []
    for index, (family, model, n, name) in enumerate(PANELS):
        y = 2-index
        points = {}
        for method in ['prior', 'ours']:
            candidates = [r for r in group(rows, family, model, method)
                          if r['singleton_true'] == r['planned']
                          and r['executed'] == r['planned']
                          and not r['node_censored'] and not r['timeout'] and not r['error']]
            if not candidates:
                continue
            best = max(candidates, key=lambda r: r['unknown_bits'])
            assert best['N'] == n
            points[method] = best['known_percentage']
            exported.append({k: best[k] for k in ['family', 'N', 'observation_model', 'method',
                                                 'unknown_bits', 'known_percentage', 'singleton_true', 'planned']})
        if len(points) == 2:
            ax.hlines(y, points['ours'], points['prior'], color='#CED5DA', lw=2.2, zorder=1)
        for method in ['prior', 'ours']:
            if method not in points:
                continue
            x = points[method]
            ax.scatter(x, y, s=62, marker=MARKERS[method], color=COLORS[method],
                       edgecolors='white', linewidth=.85, zorder=3)
            ax.annotate(f'{x:.2f}', (x, y), xytext=(0, 9), textcoords='offset points',
                        ha='center', va='bottom', fontsize=10.5, color=COLORS[method])
    ax.set_yticks([2,1,0], [f'{name} ($N={n}$)' for _,_,n,name in PANELS])
    ax.set_ylim(-.55, 2.55)
    ax.set_xlim(50,90)
    ax.set_xticks([50,60,70,80,90])
    ax.set_xlabel('Revealed bits (%)', labelpad=8)
    ax.tick_params(axis='y', length=0, pad=12)
    ax.tick_params(axis='x', length=3, color='#8A949B')
    ax.grid(axis='x', color='#E7EBEE', linewidth=.65)
    ax.set_axisbelow(True)
    ax.spines[['left','top','right']].set_visible(False)
    ax.spines['bottom'].set_color('#AAB2B8')
    ax.legend(handles=legend_handles(), loc='lower center', bbox_to_anchor=(.5,1.04),
              ncol=2, frameon=False, columnspacing=1.8, handletextpad=.35)
    save(fig, output, 'leakage')
    return exported


def recovery_plot(rows, output):
    fig = plt.figure(figsize=(9,3.6))
    gs = fig.add_gridspec(2,3, height_ratios=[.13,1], hspace=.06, wspace=.15,
                         left=.075, right=.985, top=.84, bottom=.26)
    for column, (family, model, n, name) in enumerate(PANELS):
        band = fig.add_subplot(gs[0,column])
        ax = fig.add_subplot(gs[1,column], sharex=band)
        band.set_facecolor('#F2F4F5')
        band.set_ylim(-.6,.6)
        band.set_yticks([])
        band.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        band.spines[:].set_visible(False)
        band.set_title(f'{name}, $N={n}$', fontsize=11, pad=8)
        all_x = []
        for method in ['prior', 'ours']:
            selected = group(rows, family, model, method)
            all_x.extend(r['known_percentage'] for r in selected)
            valid = [r for r in selected if r['executed'] == r['planned']
                     and not r['node_censored'] and not r['timeout'] and not r['error']]
            for a, b in zip(valid, valid[1:]):
                if b['unknown_bits'] - a['unknown_bits'] == 1:
                    ax.plot([a['known_percentage'], b['known_percentage']],
                            [100*a['singleton_true']/a['planned'], 100*b['singleton_true']/b['planned']],
                            color=COLORS[method], linewidth=1.1, zorder=2)
            ax.scatter([r['known_percentage'] for r in valid],
                       [100*r['singleton_true']/r['planned'] for r in valid],
                       color=COLORS[method], marker=MARKERS[method],
                       s=33 if method=='prior' else 23, edgecolor='white', linewidth=.4, zorder=3)
            censored = [r for r in selected if r['node_censored'] or r['timeout']
                        or r['error'] or r['executed'] < r['planned']]
            band.scatter([r['known_percentage'] for r in censored], [0]*len(censored),
                         color='#777777', marker='x', s=23, linewidths=1)
        ax.set_xlim(min(all_x)-.6, max(all_x)+.6)
        ax.set_ylim(-4,104)
        ax.set_yticks([0,50,100])
        ax.tick_params(axis='y', labelleft=column==0, length=3)
        ax.tick_params(axis='x', length=3)
        ax.grid(axis='y', color='#E5EAED', linewidth=.6)
        ax.spines[['top','right']].set_visible(False)
        ax.spines[['left','bottom']].set_color('#AAB2B8')
        ax.locator_params(axis='x', nbins=4)
    fig.text(.018,.50,'Recovery (%)',rotation=90,ha='center',va='center',fontsize=11)
    fig.supxlabel('Revealed bits (%)',y=.105,fontsize=11)
    fig.legend(handles=legend_handles(True),loc='lower center',bbox_to_anchor=(.5,.005),
               ncol=3,frameon=False,columnspacing=1.7,handletextpad=.35)
    save(fig,output,'recovery')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=ROOT/'results/summary.csv')
    parser.add_argument('--out',type=Path,default=ROOT/'fig')
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'xtick.labelsize':10,
                         'ytick.labelsize':11,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    rows=load_rows(args.data)
    leakage_plot(rows,args.out)
    recovery_plot(rows,args.out)
    print(f'Figures written to {args.out}')


if __name__=='__main__':
    main()
