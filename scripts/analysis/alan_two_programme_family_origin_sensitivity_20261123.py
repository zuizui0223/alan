#!/usr/bin/env python3
"""Independent-source audit of aquatic/terrestrial light-trap contrasts.

One 2016 Parkinson/Lawson/Tiegs mesocosm experiment and one 2022 SevenPonds
placement/colour experiment. Sensitivity analyses NEVER increment programme k.
Use source family counts (not mixed published summary totals) and whole-cell
fixed menus; no pseudocounts, excluded dark arms, or invented independence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote, urlparse
from urllib.request import urlopen, Request

ROOT = Path(__file__).resolve().parents[2]
SOURCE_BASE = ROOT / 'database/alan_experiments/parkinson2023_sources'
SOURCE_COMMIT = '6985359fd29dfe321c740691bf97ea9d564db769'
AUTHOR_SOURCE_URL = ('https://raw.githubusercontent.com/LizzParkinson/ALAN_Data/'
                     + SOURCE_COMMIT + '/Aquatic%20and%20Terrestrial%20Invertebrate%20Inputs.csv')
# All SHA1 refer to *Git blob* object IDs, not unwrapped content SHA1.
GIT_OBJECTS = {
    'parkinson_2020_authors': '4e2b7065d67e0604ca121ea26c96b5db6d64de94',
    'parkinson_2023_traps': 'c544c4af5c39677f8d5d4bacbd40ef2a518868fb',
    'parkinson_2023_families': 'd99b0e97d6608d8420342394ca5e9c65a89011c5',
}
COLOURS = ('green','purple','red','white')
EXCLUDE_BASE = ('ephydridae',)
PRIMARY_VARIANTS = (
    ('harmonized_excl_ephydridae', EXCLUDE_BASE),
    ('harmonized_excl_ephydridae_caenidae', EXCLUDE_BASE+('caenidae',)),
    ('harmonized_excl_ephydridae_leptoceridae', EXCLUDE_BASE+('leptoceridae',)),
    ('harmonized_excl_ephydridae_both_dominant', EXCLUDE_BASE+('caenidae','leptoceridae')),
)
T_975_DF1 = 12.706204736432095


def git_blob_sha(raw:bytes)->str:
    return hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()


def assert_pinned(raw:bytes, blob_id:str, source:str)->None:
    if git_blob_sha(raw)!=blob_id:
        raise ValueError('PINNED_SOURCE_BLOB_MISMATCH:'+source)


def read_author_bytes()->bytes:
    parsed=urlparse(AUTHOR_SOURCE_URL)
    if parsed.scheme!='https' or parsed.hostname!='raw.githubusercontent.com' or parsed.path.count('/')<4:
        raise ValueError('UNEXPECTED_AUTHOR_SOURCE_URL')
    # Fixed commit only. Accessing authentic public author records without
    # republishing the raw individual mesocosm rows.
    req=Request(AUTHOR_SOURCE_URL, headers={'User-Agent':'ALAN-source-audit/2026-11-23'})
    with urlopen(req,timeout=30) as resp:
        if urlparse(resp.geturl()).hostname!='raw.githubusercontent.com':
            raise ValueError('AUTHOR_SOURCE_REDIRECT_OUTSIDE_ALLOWLIST')
        raw=resp.read(100_001)
    if len(raw)>100_000:raise ValueError('AUTHOR_SOURCE_TOO_LARGE')
    assert_pinned(raw,GIT_OBJECTS['parkinson_2020_authors'],'Parkinson2020')
    return raw


def parse_csv(raw:bytes)->list[dict]:
    reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
    if not reader.fieldnames:return []
    # Original 2020 author CSV includes spaces around several column labels.
    return [{k.strip():v.strip() for k,v in row.items() if k is not None} for row in reader]


def parse_inputs(authors_2020:bytes, traps_2022:bytes, taxa_2022:bytes):
    assert_pinned(authors_2020,GIT_OBJECTS['parkinson_2020_authors'],'Parkinson2020')
    assert_pinned(traps_2022,GIT_OBJECTS['parkinson_2023_traps'],'SevenPonds2022-traps')
    assert_pinned(taxa_2022,GIT_OBJECTS['parkinson_2023_families'],'SevenPonds2022-taxa')
    p20=parse_csv(authors_2020)
    p22=parse_csv(traps_2022)
    tx=parse_csv(taxa_2022)
    if len(p20)!=20 or len(p22)!=30 or len(tx)!=193:
        raise ValueError('UNEXPECTED_SOURCE_ROW_COUNTS')
    names={k.strip():k for k in p20[0]}
    family_fields={}
    for k in p20[0]:
        m=re.fullmatch(r'# ([A-Za-z]+) \(([A-Za-z])\)', k.strip())
        if m:family_fields[m.group(1).lower()]=(k, 'aquatic' if m.group(2).isupper() else 'terrestrial')
    if len(family_fields)!=36 or not {'caenidae','leptoceridae','ephydridae'}.issubset(family_fields):
        raise ValueError('PARKINSON2020_FAMILY_FIELDS_CHANGED')
    if len({r['Mesocosm ID'] for r in p20})!=20:
        raise ValueError('DUPLICATE_2020_EXPERIMENT_UNITS')
    out=[];summary_discrepancies=[]
    for z in p20:
        fam={name:float(z[key]) for name,(key,_) in family_fields.items()}
        aq=sum(fam[name] for name,(_,origin) in family_fields.items() if origin=='aquatic')
        te=sum(fam[name] for name,(_,origin) in family_fields.items() if origin=='terrestrial')
        if float(z['# individual aquatic invertberates'])!=aq or float(z['# of individual terrestrial invertebrates'])!=te:
            summary_discrepancies.append(z['Mesocosm ID'])
        if z['Treatment'] not in ('ALAN','no ALAN'):raise ValueError('BAD_PARKINSON2020_ARM')
        out.append({'programme':'Parkinson2016','unit':z['Mesocosm ID'],'habitat':'littoral',
                    'colour':'lit' if z['Treatment']=='ALAN' else 'dark',
                    'aquatic':aq,'terrestrial':te,'family':fam})
    if sum(r['colour']=='lit' for r in out)!=10 or sum(r['colour']=='dark' for r in out)!=10:
        raise ValueError('UNBALANCED_PARKINSON2020')
    family_counts=defaultdict(dict)
    family_origin={}
    for row in tx:
        uid=row['unit_id'];name=row['family_source_spelling'];value=float(row['count'])
        if uid not in {z['unit_id'] for z in p22}:raise ValueError('UNLINKED_2022_TAXON')
        if name in family_counts[uid]:raise ValueError('DUPLICATED_FAMILY_IN_TRAP')
        family_counts[uid][name]=value
        if name in family_origin and family_origin[name]!=row['origin_as_published_figure3']:
            raise ValueError('CONFLICTING_SOURCE_ORIGIN')
        family_origin[name]=row['origin_as_published_figure3']
    assert family_origin['ephydridae']=='terrestrial'
    assert family_origin['caenidae']==family_origin['leptoceridae']=='aquatic'
    for z in p22:
        uid=z['unit_id'];fam=family_counts[uid]
        aq=sum(val for name,val in fam.items() if family_origin[name]=='aquatic')
        te=sum(val for name,val in fam.items() if family_origin[name]=='terrestrial')
        if abs(aq-float(z['aquatic_origin_capture']))>1e-10 or abs(te-float(z['terrestrial_origin_capture']))>1e-10:
            raise ValueError('ORIGIN_SUMMARY_NOT_EQUAL_FAMILY_ROWS:'+uid)
        if abs(aq+te-float(z['all_invertebrate_capture']))>1e-10:
            raise ValueError('TOTAL_NOT_EQUAL_SOURCE_ROWS:'+uid)
        out.append({'programme':'SevenPonds2022','unit':uid,'habitat':z['habitat'],
                    'colour':z['colour'],'aquatic':aq,'terrestrial':te,'family':fam})
    cells=defaultdict(int)
    for row in out:
        if row['programme']=='SevenPonds2022':cells[(row['habitat'],row['colour'])]+=1
    if len(cells)!=10 or any(n!=3 for n in cells.values()):
        raise ValueError('UNBALANCED_SEVENPONDS')
    if len({z['unit'] for z in out if z['programme']=='SevenPonds2022'})!=30:
        raise ValueError('DUPLICATE_2022_TRAPS')
    if not {'caenidae','leptoceridae'}.issubset(family_origin):raise ValueError('MISSING_2022_FAMILIES')
    return out,{'n2020':20,'n2022':30,'n2022_taxon_rows':len(tx),
                'source_defined_2020_aquatic_summary_discrepancy_units':summary_discrepancies,
                'source_2020_ephydridae_origin':'aquatic',
                'source_2022_ephydridae_origin':'terrestrial',
                'source_2020_git_blob':GIT_OBJECTS['parkinson_2020_authors'],
                'source_2022_traps_git_blob':GIT_OBJECTS['parkinson_2023_traps'],
                'source_2022_taxa_git_blob':GIT_OBJECTS['parkinson_2023_families']}


def excluded_rows(rows,excluded):
    expected={'ephydridae','caenidae','leptoceridae'}
    if not set(excluded).issubset(expected):raise ValueError('UNSAFE_EXCLUSION_CLASS')
    out=[]
    for row in rows:
        new=dict(row)
        for f in excluded:
            v=row['family'].get(f,0)
            # True source-specific conflicting Ephydridae origin is resolved
            # by excluding that family from both programmes.
            key='terrestrial' if f=='ephydridae' and row['programme']=='SevenPonds2022' else 'aquatic'
            new[key]-=v
        if any(new[k]<0 for k in ('aquatic','terrestrial')):
            raise ValueError('NEGATIVE_AFTER_SOURCE_EXCLUSION')
        out.append(new)
    return out


def avg(x):
    if not x:raise ValueError('EMPTY_MEAN')
    return sum(x)/len(x)


def sample_cov(x,y):
    if len(x)!=len(y) or len(x)<2:raise ValueError('NO_INDEPENDENT_UNITS')
    ma,mb=avg(x),avg(y)
    return sum((a-ma)*(b-mb) for a,b in zip(x,y))/(len(x)-1)


def one_programme(rows,programme,omit_colour=None):
    selected=[r for r in rows if r['programme']==programme]
    if not selected:return None
    arms={};trace={}
    for is_lit,label in [(True,'lit'),(False,'dark')]:
        arm=[r for r in selected if (r['colour']!='dark')==is_lit and not (is_lit and r['colour']==omit_colour)]
        cell=defaultdict(list)
        for r in arm:cell[(r['habitat'],r['colour'])].append(r)
        groups=list(cell.values());n_groups=len(groups)
        if not groups:return None
        ma=avg([avg([r['aquatic'] for r in g]) for g in groups])
        mt=avg([avg([r['terrestrial'] for r in g]) for g in groups])
        if ma<=0 or mt<=0:raise ValueError('NONPOSITIVE_ARM_MEAN_STOP')
        v=0
        for g in groups:
            aq=[r['aquatic'] for r in g];te=[r['terrestrial'] for r in g];n=len(g)
            v+=(sample_cov(aq,aq)/ma**2+sample_cov(te,te)/mt**2-2*sample_cov(aq,te)/(ma*mt))/(n_groups*n_groups*n)
        arms[label]={'aquatic_mean':ma,'terrestrial_mean':mt,'variance':v,
                     'physical_light_units':sum(len(g) for g in groups),'treatment_cells':n_groups}
    logrr=math.log(arms['lit']['aquatic_mean']/arms['lit']['terrestrial_mean'])-math.log(arms['dark']['aquatic_mean']/arms['dark']['terrestrial_mean'])
    return {'log_ratio':logrr,'ratio':math.exp(logrr),'sampling_variance':arms['lit']['variance']+arms['dark']['variance'],
            'arms':arms}


def two_programme_mkh(a,b):
    ys=(a['log_ratio'],b['log_ratio']);vs=(a['sampling_variance'],b['sampling_variance'])
    if not all(v>0 for v in vs):raise ValueError('NONPOSITIVE_ESTIMATED_VARIANCE')
    tau2=max(0,((ys[0]-ys[1])**2-sum(vs))/2)
    weights=[1/(v+tau2) for v in vs]
    mu=sum(w*y for w,y in zip(weights,ys))/sum(weights)
    q=sum(w*(y-mu)**2 for w,y in zip(weights,ys))
    variance=max(1,q)/sum(weights)
    radius=T_975_DF1*math.sqrt(variance)
    return {'ratio':math.exp(mu),'log_ratio':mu,'lower':math.exp(mu-radius),'upper':math.exp(mu+radius),
            'tau2_REML_two_programme_closed_form':tau2,'Q':q,'df':1,'w_2020':weights[0]/sum(weights),
            'mKH_sampling_variance':variance,'programmes':2}


def source_family_counts(rows,family):
    output={}
    for programme in ('Parkinson2016','SevenPonds2022'):
        cohort=[r for r in rows if r['programme']==programme]
        lit=[r for r in cohort if r['colour']!='dark'];dark=[r for r in cohort if r['colour']=='dark']
        output[programme]={'lit':sum(r['family'].get(family,0) for r in lit),
                           'dark':sum(r['family'].get(family,0) for r in dark)}
    return output


def analyse(rows,source_audit):
    variants=[];benchmark=None
    for label,removed in PRIMARY_VARIANTS:
        adjusted=excluded_rows(rows,removed)
        a=one_programme(adjusted,'Parkinson2016');b=one_programme(adjusted,'SevenPonds2022')
        fit=two_programme_mkh(a,b)
        variants.append({'variant':label,'excluded_families':list(removed),
                         'Parkinson2016':a,'SevenPonds2022':b,'pooled':fit})
        if removed==EXCLUDE_BASE:benchmark=(a,b,fit)
    if benchmark is None:raise AssertionError('baseline missing')
    b20,b22,bpooled=benchmark
    if abs(b20['ratio']-5.12503359312013)>1e-10 or abs(b22['ratio']-12.175906030654636)>1e-10 or abs(bpooled['ratio']-6.88025420158834)>1e-10:
        raise ValueError('EXISTING_TWO_PROGRAMME_META_DID_NOT_REPLAY')
    spectral=[]
    normalized=excluded_rows(rows,EXCLUDE_BASE)
    for colour in COLOURS:
        b=one_programme(normalized,'SevenPonds2022',omit_colour=colour)
        spectral.append({'omitted_colour':colour,'Parkinson2016':b20,
                         'SevenPonds2022':b,'pooled_menu_sensitivity':two_programme_mkh(b20,b),
                         'status':'unstandardized_wavelength_menu_perturbation_not_new_study'})
    audit={'canonical_effect_catalogue_altered':False,'canonical_experimental_k_added':0,
           'programme_count':2,'source_pin':source_audit,'family_counts':{
               'caenidae':source_family_counts(rows,'caenidae'),
               'leptoceridae':source_family_counts(rows,'leptoceridae'),
               'ephydridae':source_family_counts(rows,'ephydridae')},
           'distinct_2020_2022_experiments':'two Michigan freshwater-edge programmes, differing year and equipment',
           'endpoint':'ratio of aquatic-origin to terrestrial-origin *trap catch* responses; not emergence, arthropod abundance, predator attack, pollination or plant fitness',
           'inference':'post hoc dominant-family and colour-menu sensitivity, k=2; confidence bounds modified Knapp-Hartung df=1 and not definitive',
           'limitations':['source 2020 summed taxa mismatch author summary in four mesocosms',
                          'Ephydridae different published origin classes, excluded in harmonized analysis',
                          'colour spectra/illuminance and habitat positions differ between programmes',
                          'ponds are both in Michigan, no independent wet-vs-dry plant habitat comparison',
                          'one experiment 2022 only single pond one night; multiple traps can attract one moving insect population',
                          'high dependence of aquatic capture response on Caenidae and Leptoceridae']}
    return {'source':audit,'variants':variants,'leave_one_colour_out':spectral}


def csv_rows(doc):
    out=[]
    for row in doc['variants']:
        a,b,m=row['Parkinson2016'],row['SevenPonds2022'],row['pooled']
        out.append({'variant':row['variant'],'removed_families':'+'.join(row['excluded_families']),
                    'RR_2016':a['ratio'],'V_2016':a['sampling_variance'],
                    'RR_2022':b['ratio'],'V_2022':b['sampling_variance'],
                    'pooled_RR':m['ratio'],'lower_mKH':m['lower'],'upper_mKH':m['upper'],
                    'tau2':m['tau2_REML_two_programme_closed_form'],
                    'k_independent_programmes':2,'status':'exploratory_not_added_to_canonical_catalogue'})
    for row in doc['leave_one_colour_out']:
        a,b,m=row['Parkinson2016'],row['SevenPonds2022'],row['pooled_menu_sensitivity']
        out.append({'variant':'2022_exclude_'+row['omitted_colour']+'_colour_cell',
                    'removed_families':'ephydridae','RR_2016':a['ratio'],'V_2016':a['sampling_variance'],
                    'RR_2022':b['ratio'],'V_2022':b['sampling_variance'],
                    'pooled_RR':m['ratio'],'lower_mKH':m['lower'],'upper_mKH':m['upper'],
                    'tau2':m['tau2_REML_two_programme_closed_form'],
                    'k_independent_programmes':2,'status':'different_colour_menu_exploratory_not_comparable_common_wavelength'})
    return out


def save(doc,dir:Path):
    dir.mkdir(parents=True,exist_ok=True)
    (dir/'source_harmonization_sensitivity.json').write_text(json.dumps(doc,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    rows=csv_rows(doc)
    with (dir/'family_and_colour_sensitivity.csv').open('w',encoding='utf-8',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',default=str(ROOT/'database/alan_experiments/arthropod_expansion_20261122/two_programme_taxon_sensitivity'))
    args=parser.parse_args()
    authors=read_author_bytes()
    trap=(SOURCE_BASE/'trap_units.csv').read_bytes()
    taxa=(SOURCE_BASE/'taxon_counts.csv').read_bytes()
    rows,audit=parse_inputs(authors,trap,taxa)
    doc=analyse(rows,audit);save(doc,Path(args.output_dir))
    print('VERIFIED_TWO_PROGRAMME_TAXON_SENSITIVITY',len(doc['variants']),len(doc['leave_one_colour_out']),len(rows))
    for v in doc['variants']:
        print(v['variant'],f"2016={v['Parkinson2016']['ratio']:.4f}",f"2022={v['SevenPonds2022']['ratio']:.4f}",f"pooled={v['pooled']['ratio']:.4f}")


if __name__=='__main__':main()
