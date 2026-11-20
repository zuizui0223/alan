"""Independent re-extraction: standard library only, no project helpers imported."""
import csv
import hashlib
import json
import math
import statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = REPO.parent / 'gbif_annual/alan_experiment_sources/parkinson2020'
DB = HERE.parent

def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def write(name, rows):
    with (HERE/name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

pan = read(RAW/'Aquatic and Terrestrial Invertebrate Inputs.csv')
web = read(RAW/'Spider and Prey Insect Abundance.csv')
assert len(pan) == 20 and len(web) == 60
assert {r['Treatment'] for r in pan+web} == {'ALAN', 'no ALAN'}
assert len({r['Mesocosm ID'] for r in pan}) == 20
ids = ['P2020_'+s for s in ['aquatic_pan','terrestrial_pan','spider_final','prey_final','spider_mean_3dates','prey_in_web_mean_3dates']]
units = []
for p in sorted(pan, key=lambda r:int(r['Mesocosm ID'])):
    rows = [r for r in web if r['Mesocosm ID'] == p['Mesocosm ID']]
    assert len(rows) == 3 and {r['Sampling Date'] for r in rows} == {'7/8/2016','7/23/2016','8/4/2016'}
    assert all(r['Treatment'] == p['Treatment'] for r in rows)
    last = next(r for r in rows if r['Sampling Date'] == '8/4/2016')
    vals = [float(p['# individual aquatic invertberates']),float(p['# of individual terrestrial invertebrates']),
            float(last['# Individual tetragnathids']),float(last['# individual insects in webs']),
            st.mean(float(r['# Individual tetragnathids']) for r in rows),st.mean(float(r['# individual insects in webs']) for r in rows)]
    units.append({'mesocosm_id':int(p['Mesocosm ID']), 'treatment':p['Treatment'], **dict(zip(ids, vals))})
write('unit_values.csv', units)
groups = {g:[u for u in units if u['treatment']==g] for g in ['ALAN','no ALAN']}
assert all(len(v)==10 for v in groups.values())
existing = {r['effect_id']:r for r in read(DB/'effects_computed.csv')}
effects, comparisons = [], []
for k in ids:
    x,y = ([r[k] for r in groups[g]] for g in ['ALAN','no ALAN'])
    mx,my=st.mean(x),st.mean(y)
    sx,sy=st.stdev(x),st.stdev(y)
    out={'effect_id':k,'mean_lit':mx,'mean_dark':my,'dispersion_lit':sx,'dispersion_dark':sy,'n_lit':len(x),'n_dark':len(y),
         'lnRR':math.log(mx/my),'sampling_variance':sx*sx/(len(x)*mx*mx)+sy*sy/(len(y)*my*my)}
    effects.append(out)
    for field,value in out.items():
        if field=='effect_id': continue
        old=float(existing[k][field])
        comparisons.append({'effect_id':k,'field':field,'independent':value,'existing':old,'difference':value-old,'pass':math.isclose(value,old,rel_tol=1e-11,abs_tol=1e-12)})
write('effects_reextracted.csv',effects)
write('effect_comparison.csv',comparisons)
previous={(r['effect_id_1'],r['effect_id_2']):float(r['sampling_covariance']) for r in read(DB/'effect_covariance.csv')}
covrows=[]
for a in ids:
    for b in ids:
        cov=0.0
        for g,rs in groups.items():
            xa=[r[a] for r in rs];xb=[r[b] for r in rs]
            ma,mb=st.mean(xa),st.mean(xb);n=len(rs)
            sample_cov=sum((u-ma)*(v-mb) for u,v in zip(xa,xb))/(n-1)
            cov+=sample_cov/(n*ma*mb)
        old=previous[(a,b)]
        covrows.append({'effect_id_1':a,'effect_id_2':b,'independent_covariance':cov,'existing_covariance':old,'difference':cov-old,'pass':math.isclose(cov,old,rel_tol=1e-11,abs_tol=1e-12)})
write('covariance_comparison.csv',covrows)
# Preserve classification discrepancies rather than silently recalculating origin totals.
anomalies=[]
headers=list(pan[0]); split=headers.index('# Cosmopterygidae (a)')
for p in pan:
    aq=sum(int(p[k]) for k in headers[6:split]);te=sum(int(p[k]) for k in headers[split:])
    if aq!=int(p[headers[4]]) or te!=int(p[headers[5]]):
        anomalies.append({'mesocosm_id':p['Mesocosm ID'],'aquatic_as_reported':int(p[headers[4]]),'aquatic_family_sum':aq,'terrestrial_as_reported':int(p[headers[5]]),'terrestrial_family_sum':te})
summary={'all_effect_fields_pass':all(r['pass'] for r in comparisons),'all_36_covariances_pass':all(r['pass'] for r in covrows),
         'max_effect_absolute_difference':max(abs(r['difference']) for r in comparisons),'max_covariance_absolute_difference':max(abs(r['difference']) for r in covrows),
         'origin_family_sum_discrepancies':anomalies,'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in RAW.glob('*.csv')},
         'verification_type':'independent_agent_computational_reextraction; not human review'}
(HERE/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
