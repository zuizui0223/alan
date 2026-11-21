"""Render the fixed-species composition result; no pooled-study diamond."""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def main(source, out):
    out.parent.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(source/'composition_and_trophic_contrasts.csv')
    species = pd.read_csv(source/'species_marginal_effects.csv')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    colours = {'L2': '#B56C13', 'L3': '#2563A6'}
    offsets = {'L2': .13, 'L3': -.13}
    fig = plt.figure(figsize=(15.5, 7.2), facecolor='white')
    gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.05, 1.23],
                          left=.16, right=.985, bottom=.29, top=.77, wspace=.67)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

    def draw(ax, subset, y, intervals='nominal', marker='o'):
        for r in subset.itertuples():
            yy = y+offsets[r.light]
            lo = r.lower_95_t5 if intervals == 'nominal' else r.lower_95_family_t5
            hi = r.upper_95_t5 if intervals == 'nominal' else r.upper_95_family_t5
            ax.errorbar(r.ratio, yy, xerr=[[r.ratio-lo], [hi-r.ratio]],
                fmt=marker, color=colours[r.light], markersize=6, capsize=3,
                linewidth=1.5, markeredgecolor='white', markeredgewidth=.5, zorder=3)

    names = ['Agrostis tenuis', 'Anthoxanthum odoratum', 'Hypochaeris radicata', 'Lotus pedunculatus']
    traits = pd.read_csv(source/'penryn_four_species_PLANTATT_verified.tsv', sep='\t').set_index('source_species')
    labels = [f'{name}\nF = {int(traits.loc[name, "F_PLANTATT"])}' for name in names]
    for i, name in enumerate(names):
        draw(axes[0], species[species.species.eq(name)], 3-i)
    axes[0].set(yticks=range(4), yticklabels=labels[::-1], ylim=(-.55, 3.55),
                xlim=(.32, 3.5), xlabel='Biomass response ratio\nALAN / dark')
    axes[0].set_title('A  Species responses', loc='left', fontweight='bold', pad=24)
    axes[0].text(0, 1.025, 'Equal weight across T1, T2 and T3', transform=axes[0].transAxes,
                 fontsize=9, color='#555555')

    for i, trophic in enumerate(['T1', 'T2', 'T3']):
        draw(axes[1], data[data.trophic.eq(trophic)], 3-i)
    draw(axes[1], data[data.kind.eq('marginal')], 0, 'family', marker='s')
    axes[1].axhline(.55, color='#dddddd', linewidth=.8)
    axes[1].set(yticks=range(4), yticklabels=['Average\n(2-colour intervals)', 'Plants + herbivores\n+ predators (T3)',
        'Plants + herbivores\n(T2)', 'Plants only\n(T1)'], ylim=(-.55,3.55),
        xlim=(.075, 3), xlabel='Relative response ratio\nLotus / three comparator responses')
    axes[1].set_title('B  Relative response', loc='left', fontweight='bold', pad=24)
    axes[1].text(0, 1.025, 'Same-unit species covariance retained', transform=axes[1].transAxes,
                 fontsize=9, color='#555555')

    draw(axes[2], data[data.trophic.eq('T2_minus_T1')], 1, 'family')
    draw(axes[2], data[data.trophic.eq('T3_minus_T2')], 0, 'family')
    axes[2].set(yticks=[0,1], yticklabels=['Predator addition\nT3 vs T2', 'Herbivore package\nT2 vs T1'],
        ylim=(-.55,1.55), xlim=(.075,20), xlabel='Modification of relative response\nRatio of the two light-response contrasts')
    axes[2].set_title('C  Trophic modification', loc='left', fontweight='bold', pad=24)
    axes[2].text(0, 1.025, '4-comparison intervals', transform=axes[2].transAxes,
                 fontsize=9, color='#555555')

    for i, ax in enumerate(axes):
        ax.set_xscale('log')
        ax.axvline(1, color='#7b8794', linestyle='--', linewidth=1, zorder=1)
        ax.grid(axis='x', which='major', color='#e7ebee', linewidth=.65)
        ticks = [.5,1,2] if i==0 else ([.1,.25,.5,1,2] if i==1 else [.1,.25,.5,1,2,5,10])
        ax.set_xticks(ticks, [str(x) for x in ticks])
        ax.tick_params(axis='y', length=0, pad=8)
        ax.spines[['top','right','left']].set_visible(False)
        ax.spines['bottom'].set_color('#b0b8c0')
        ax.minorticks_off()
    fig.suptitle('Penryn: Lotus biomass responses relative to three other species\n'
                 'within the same lighting and trophic experiment', x=.05, y=.97,
                 ha='left', fontsize=18, fontweight='bold', color='#182c3b')
    fig.legend(handles=[Line2D([0],[0],color=colours['L2'],marker='o',label='Amber LED'),
                        Line2D([0],[0],color=colours['L3'],marker='o',label='White LED')],
               loc='upper right', bbox_to_anchor=(.985,.905), ncol=2, frameon=False)
    fig.text(.05,.14,'One experiment | 54 independently illuminated or dark mesocosms | 4 complete species | 6 units per treatment cell',
             fontsize=10, color='#34495e')
    fig.text(.05,.09,'Intervals: approximate t5; panel A and B cells nominal 95%; B averages Bonferroni for 2 colours; C Bonferroni for 4 modifiers.',
             fontsize=9, color='#555555')
    fig.text(.05,.052,'Post hoc subset. Only Lotus has F > 7; species identity and moisture affinity cannot be separated. T1 includes insecticide/molluscicide.',
             fontsize=9, color='#555555')
    for extension in ['png','pdf','svg']:
        fig.savefig(out.with_suffix('.'+extension), dpi=220,
                    metadata={'Creator':'alan reproducible analysis'} if extension in ['pdf','svg'] else None)
    plt.close(fig)
    svg_path = out.with_suffix('.svg')
    svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines()) + '\n')
    manifest = dict(inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [
        source/'composition_and_trophic_contrasts.csv', source/'species_marginal_effects.csv',
        source/'penryn_four_species_PLANTATT_verified.tsv']},
        displayed_estimates=len(species)+len(data), covariance_aware=True, pooled_study_estimate=False,
        figures=[out.with_suffix('.'+s).name for s in ['png','pdf','svg']])
    out.with_suffix('.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args(); main(args.source,args.out)
