"""Plot reported S4b contrasts, with approximate t intervals; no meta pooling."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import t

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parents[2] / 'manuscript/biotic_synthesis/figures'
OUT.mkdir(parents=True, exist_ok=True)
source = ROOT / 'table_S4b.csv'
with source.open(encoding='utf-8', newline='') as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 12
derived = []
species = response = None
for row in rows:
    species = row['Species'] or species
    response = row['Response'].strip() or response
    assert species and response
    assert row['Contrast'] == 'Unlit - Lit'
    estimate, se, df = -float(row['Estimate']), float(row['SE']), float(row['df'])
    margin = t.ppf(.975, df) * se
    derived.append(dict(species=species, response=response,
        access=row['Pollination treatment'], direction='Lit minus Unlit',
        estimate=estimate, SE=se, df=df, ci_lower=estimate-margin, ci_upper=estimate+margin,
        reported_p=row['p value'], source_table='S4b',
        interval_method='approximate marginal t interval from rounded reported SE and df; no multiplicity adjustment',
        scale='reported model scale; main Methods not verified',
        pooling_status='not_pooled_independence_and_scale_unverified'))
with (ROOT / 'reported_contrasts_lit_minus_unlit.csv').open('w', encoding='utf-8', newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(derived[0]));writer.writeheader();writer.writerows(derived)
fig, axes=plt.subplots(1,2,figsize=(11,5.6))
for ax,response in zip(axes,['Pollen number','Seed set']):
    selected=[r for r in derived if r['response']==response]
    assert len(selected)==6
    for y,r in enumerate(selected):
        ax.errorbar(r['estimate'],y,xerr=[[r['estimate']-r['ci_lower']],[r['ci_upper']-r['estimate']]],
            fmt='o',color='#215a8a' if r['species']=='P. kansuensis' else '#ad6429',capsize=3)
    ax.axvline(0,color='0.5',linestyle='--',linewidth=1)
    ax.set_yticks(range(6),[r['species']+' | '+r['access'] for r in selected],fontsize=9)
    ax.invert_yaxis();ax.set_title(response)
    ax.set_xlabel('Lit minus Unlit (reported model scale)')
    ax.spines[['top','right']].set_visible(False)
fig.suptitle('Reported lighting contrasts by pollinator-access condition',fontsize=13)
fig.text(.5,.03,'One report, dependent contrasts; approximate marginal 95% t intervals.\nMain Methods and seed-set scale unverified; not a meta-analysis or percent change.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.11,1,.94))
for ext in ['png','pdf','svg']:
    fig.savefig(OUT / ('lu_jpe_access_contrasts.'+ext),dpi=180)
svg = OUT / 'lu_jpe_access_contrasts.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n', encoding='utf-8')
audit=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    source_rows=12, derived_rows=len(derived), sign_change='Estimate negated; SE/df unchanged',
    label_handling='Blank species/response cells forward-filled within S4b only; source table CSV unchanged',
    intervals='t(0.975, reported df) times reported SE; approximate, marginal, not adjusted',
    no_meta_analysis=True, no_between_species_test=True, covariance_not_imputed=True)
(ROOT / 'reported_contrast_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
