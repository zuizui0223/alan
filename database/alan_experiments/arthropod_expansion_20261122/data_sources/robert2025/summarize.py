"""Reproduce paired descriptive captures from Robert et al. Figshare worksheet Insects_Order."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
from statistics import mean, stdev, covariance
import openpyxl

root=Path(__file__).resolve().parent
file=root/'source.xlsx'
meta=json.loads((root/'metadata.json').read_bytes())
assert hashlib.md5(file.read_bytes()).hexdigest()==meta['files'][0]['computed_md5']
book=openpyxl.load_workbook(file, read_only=True, data_only=True)
all_rows=list(book['Insects_Order'].values)
rows=[row for row in all_rows[1:] if row[0] is not None]
assert len(rows)==96 and len(all_rows)==98
cells=defaultdict(list)
for row in rows:
    datestr,site,arm,time,abundance=row[:5]
    assert sum(row[6:14])==abundance
    cells[site,arm,time].append(abundance)
assert len(cells)==32 and all(len(v)==3 for v in cells.values())
assert sum(row[4] for row in rows)==all_rows[-1][4]==375
sites=sorted({site for site,arm,time in cells})
def save(name, records):
    with (root/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
site_rows=[]; contrasts=[];values={}
for site in sites:
    for arm in ('C','A','FW','CW','all_lit_site_mean'):
        treatments=('A','FW','CW') if arm=='all_lit_site_mean' else (arm,)
        before=mean(mean(cells[site,trt,'Before']) for trt in treatments)
        during=mean(mean(cells[site,trt,'During']) for trt in treatments)
        d_before=mean(cells[site,'C','Before']); d_during=mean(cells[site,'C','During'])
        diff=(during-before)-(d_during-d_before)
        site_rows.append(dict(site=site,arm=arm,n_physical_plots=len(treatments),n_nights_per_phase=3,before_captures_per_night=before,
                              during_captures_per_night=during,own_change=during-before,same_site_control_change=d_during-d_before,
                              paired_difference_in_changes=diff))
        values[site,arm]=diff
for arm in ('A','FW','CW','all_lit_site_mean'):
    data=[values[site,arm] for site in sites]
    considered=[r for r in site_rows if r['arm']==arm]
    contrasts.append(dict(arm=arm,comparator='C',n_paired_reserves=len(sites),before_arm_mean=mean(r['before_captures_per_night'] for r in considered),
                          during_arm_mean=mean(r['during_captures_per_night'] for r in considered),paired_difference_in_changes_mean=mean(data),
                          paired_difference_in_changes_sample_sd=stdev(data),paired_difference_in_changes_se=stdev(data)/(len(sites)**0.5),
                          status='descriptive_same_family_shared_controls_no_catalogue_adoption'))
save('site_capture_summaries.csv',site_rows)
save('paired_capture_checks.csv',contrasts)
cov=[]
for a,b in product(('A','FW','CW'),repeat=2):
    cov.append(dict(treatment_a=a,treatment_b=b,n_reserves=len(sites),sample_covariance_of_paired_estimates=covariance([values[s,a] for s in sites],[values[s,b] for s in sites])/len(sites)))
save('paired_capture_estimate_covariance.csv',cov)
print(contrasts)
