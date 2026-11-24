"""Rebuild week-to-plot, then plot-to-arm, descriptive capture summaries from EnviDat."""
import csv
import json
from collections import Counter, defaultdict
from datetime import date, timedelta
from pathlib import Path
from statistics import mean, stdev

root = Path(__file__).resolve().parent
def read(name):
    with (root/name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))
def write(name, rows):
    with (root/name).open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
def export_summary(acc, endpoint, period, status, site, arm, subset):
    values=[row['mean_captures_per_sample_week'] for row in subset]
    acc.append(dict(endpoint=endpoint,period=period,scope=status,site=site,arm=arm,
                    n_plots=len(values),n_plot_week_samples=sum(r['n_plot_week_samples'] for r in subset),
                    unweighted_mean_of_plot_means=mean(values),sample_sd_across_plot_means=stdev(values)))

taxon_fields={'flight_interception':'Nematocera Brachycera Lepidoptera Coleoptera Hymenoptera Hemiptera Neuropteroidea Trichoptera',
              'ground_pitfall':'Collembola Acarina Formicidae Arachnida Coleo_Pred Coleo_Rest Myriapoda Hemiptera'}
output, stats, issues, overlap, coords = [], [], [], {}, []
for endpoint, code in [('flight_interception','pt'),('ground_pitfall','trf')]:
    raw=read(f'566/alanex_{code}_all_clean_control.csv')
    used=[]
    for i,row in enumerate(raw,2):
        if (row['temp'] not in ('2000','3000','4000','control') or row['level'] not in ('50','100','control') or row['issues']!='0'):
            continue
        assert sum(int(row[k]) for k in taxon_fields[endpoint].split())==int(row['SumAll'])
        used.append(row)
        if code=='pt':
            dat=date(1899,12,30)+timedelta(days=int(row['CollectionDate']))
            tags=[]
            if str(dat.year)!=row['year']:tags.append('collection_date_year_mismatch')
            if row['plot']=='L15':tags.append('extra_plot_in_retained_arm')
            if row['plot']=='L12' and row['year']=='2022' and row['week_ID']=='3':tags.append('one_lamp_two_control_states_same_week')
            if tags:
                issues.append(dict(source=f'alanex_{code}_all_clean_control.csv',source_row=i,site=row['site'],
                                   plot=row['plot'],year=row['year'],week_ID=row['week_ID'],collection_date_decoded=dat.isoformat(),
                                   control=row['control'],SumAll=row['SumAll'],flags='|'.join(tags)))
    if code=='trf':
        older=read('496/alanex_trf_all_clean.csv')
        idx=defaultdict(list)
        for r in older:idx[(r['Site'],r['Plot'],r['CollectionDate'],r['Position'])].append(r)
        cols='Collembola Acarina Formicidae Carabidae Araneae Opiliones Chilopoda'.split()
        exact=0
        mismatches=[]
        for r in used:
            key=(r['site'],r['plot'],r['CollectionDate'],r['position'])
            candidates=idx[key]
            if len(candidates)==1 and all(r[c]==candidates[0][c] for c in cols):exact+=1
            else:mismatches.append(list(key))
        overlap={'comparison':'same site, plot, collection date and pitfall position; agreement for each named elementary taxon',
                 'taxa':cols,'eligible_566_ground_rows':len(used),'exact_matches':exact,'unmatched':mismatches,
                 'warning':'The archives define SumAll differently; only 2/306 shared center samples have equal SumAll. Do not use both reports as independent data.'}
    byplot=defaultdict(list)
    for r in used:
        arm='dark' if r['control']=='y' else 'lit'
        byplot[(r['site'],r['plot'],arm,r['year'])].append(r)
    per_year=[]
    for (site,plot,arm,year),records in sorted(byplot.items()):
        week_ids=sorted(int(r['week_ID']) for r in records)
        datum=dict(endpoint=endpoint,period=year,site=site,plot=plot,arm=arm,
                   n_plot_week_samples=len(records),sampled_week_ids=','.join(map(str,week_ids)),
                   mean_captures_per_sample_week=mean(int(r['SumAll']) for r in records))
        per_year.append(datum)
    output.extend(per_year)
    periods=sorted({r['period'] for r in per_year})
    if code=='trf':
        by_id=defaultdict(list)
        for r in per_year:by_id[(r['site'],r['plot'],r['arm'])].append(r)
        aggregate=[]
        for (site,plot,arm),data in sorted(by_id.items()):
            assert {r['period'] for r in data}=={'2021','2022'}
            aggregate.append(dict(endpoint=endpoint,period='2021-2022',site=site,plot=plot,arm=arm,
                n_plot_week_samples=sum(r['n_plot_week_samples'] for r in data),sampled_week_ids=';'.join(r['period']+':'+r['sampled_week_ids'] for r in data),
                mean_captures_per_sample_week=mean(r['mean_captures_per_sample_week'] for r in data)))
        output.extend(aggregate)
        periods+=['2021-2022']
    for year in periods:
        if code=='pt' and year=='2022':continue
        subset=[r for r in output if r['endpoint']==endpoint and r['period']==year]
        status='descriptive_observed_weeks' if not (code=='pt' and year=='2023') else 'descriptive_nominal_year_date_conflict'
        for site in ['ALL']+sorted({r['site'] for r in subset}):
            for arm in ['lit','dark']:
                cell=[r for r in subset if r['arm']==arm and site in ('ALL',r['site'])]
                assert len(cell)==({'ALL':{'lit':36,'dark':6}}.get(site,{'lit':12,'dark':2})[arm])
                export_summary(stats,endpoint,year,status,site,arm,cell)
coord_index=defaultdict(set)
for r in read('566/alanex_moransi.csv'):
    if r['site'] and r['plot'] and r['long'] and r['lat']:
        coord_index[(r['site'],r['plot'])].add((r['long'],r['lat']))
for (site,plot),locations in sorted(coord_index.items()):
    for east,north in sorted(locations):
        coords.append(dict(site=site,plot=plot,source_field_long=east,source_field_lat=north,interpretation='projected Swiss grid-looking coordinates; confirm CRS before converting'))
write('plot_summaries.csv',output)
write('arm_summaries.csv',stats)
write('flagged_flight_records.csv',issues)
write('native_coordinate_values.csv',coords)
(root/'ground_family_overlap.json').write_text(json.dumps(overlap,indent=2)+'\n')
for r in stats:
    if r['site']=='ALL':print(r['endpoint'],r['period'],r['arm'],r['n_plots'], r['n_plot_week_samples'],round(r['unweighted_mean_of_plot_means'],3),round(r['sample_sd_across_plot_means'],3))
print('flight flags',Counter(row['flags'] for row in issues))
print('overlap',overlap)
