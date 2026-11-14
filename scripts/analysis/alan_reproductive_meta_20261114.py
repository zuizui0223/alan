#!/usr/bin/env python3
"""Endpoint-specific secondary synthesis of published Guenat 2023 effects.

Uses archived SMD/SE, not individual plants and not a reproduction of the
original authors' rma.mv model. First average effects within each paper and
endpoint. Sampling covariance is unknown: use explicit rho=0,.5,.9,1 scenarios.
Then study-level REML and modified Knapp-Hartung t intervals. k<3 not pooled.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import json,math
from pathlib import Path

ENDPOINTS=('seed set','fruit set','seed number')

def study_effects(rows,rho):
    if not 0<=rho<=1:raise ValueError('Correlation outside sensitivity range')
    endpoints={r['Measurement'].strip() for r in rows}
    if len(endpoints)>1:raise ValueError('Different reproductive endpoints cannot be pooled')
    groups=defaultdict(list)
    for r in rows:
        y=float(r['SMD']);se=float(r['SE'])
        if not math.isfinite(y) or not math.isfinite(se) or se<=0:
            raise ValueError('Invalid archived effect or standard error')
        groups[r['Study'].strip()].append((y,se,r))
    out=[]
    for study,items in sorted(groups.items()):
        n=len(items);ses=[s for _,s,_ in items]
        v=((1-rho)*sum(s*s for s in ses)+rho*sum(ses)**2)/n**2
        out.append({'study':study,'yi':sum(y for y,_,_ in items)/n,'vi':v,
                    'n_source_effects':n,'doi':sorted({x['DOI'] for _,_,x in items}),
                    'rho_assumed_not_measured':rho})
    return out

def random_effects(studies):
    k=len(studies)
    if k<3:return {'status':'TOO_FEW_STUDIES_FOR_POOLED_INFERENCE','k_papers':k}
    import numpy as np
    from scipy.optimize import minimize_scalar
    from scipy.stats import t
    y=np.array([s['yi'] for s in studies],float)
    v=np.array([s['vi'] for s in studies],float)
    if not np.all(np.isfinite(y)) or not np.all(v>0):raise ValueError('Invalid study estimates')
    def objective(tau):
        w=1/(v+tau);mu=float(sum(w*y)/sum(w))
        return .5*(sum(np.log(v+tau))+math.log(sum(w))+sum(w*(y-mu)**2))
    upper=max(1.,float(np.var(y,ddof=1)*20+max(v)*20))
    fit=minimize_scalar(objective,bounds=(0,upper),method='bounded',options={'xatol':1e-12})
    if not fit.success:raise RuntimeError('REML optimization failed')
    tau=0. if objective(0)<=fit.fun else float(fit.x)
    if tau>.99*upper:raise RuntimeError('REML boundary requires a wider search')
    w=1/(v+tau);mu=float(sum(w*y)/sum(w))
    q=float(sum(w*(y-mu)**2)/(k-1))
    se=math.sqrt(max(1.,q)/sum(w));crit=float(t.ppf(.975,k-1))
    return {'status':'EXPLORATORY_STUDY_LEVEL_REML_MODIFIED_KH',
            'k_papers':k,'pooled_SMD':mu,'tau2':tau,'modified_KH_SE':se,
            'ci95':[mu-crit*se,mu+crit*se],
            'nominal_two_sided_p':float(2*t.sf(abs(mu/se),k-1)),
            'df':k-1,'HK_variance_multiplier':max(1.,q),
            'method':'paper-level REML; modified Knapp-Hartung variance floor 1; t(k-1)',
            'unmodeled_design_and_taxonomic_dependence_remain':True}

def summarize(rows):
    out={'source_effect_rows':len(rows),'papers':len({r['Study'] for r in rows}),
         'study_counts':dict(sorted(Counter(r['Study'] for r in rows).items())),
         'measurement_counts':dict(Counter(r['Measurement'] for r in rows)),
         'field_lab_counts':dict(Counter(r['Field.lab'] for r in rows)),
         'species':sorted({r['Species'] for r in rows})}
    out['by_endpoint']={}
    for ep in ENDPOINTS:
        selected=[r for r in rows if r['Measurement'].strip()==ep]
        scenarios={}
        for rho in (0.,.5,.9,1.):
            papers=study_effects(selected,rho)
            scenarios[str(rho)]={'paper_estimates':papers,'model':random_effects(papers)}
        primary=scenarios['0.5']['paper_estimates']
        deleted=[]
        if len(primary)>=4:
            for leave in primary:
                z=random_effects([s for s in primary if s['study']!=leave['study']])
                deleted.append({'omitted':leave['study'],**z})
        out['by_endpoint'][ep]={'n_effects':len(selected),'n_papers':len(primary),
            'correlation_scenarios':scenarios,'leave_one_paper_out_rho_0p5':deleted}
    return out

def main(out):
    from scripts.analysis.alan_pollution_meta_source_20261114 import load_source
    fields,rows,receipt=load_source()
    if len(rows)!=486:raise ValueError('Original meta source size changed')
    if len({r['NR'] for r in rows})!=486:raise ValueError('Duplicate source effect keys')
    alan=[r for r in rows if r['Pollutant'].strip()=='ALAN']
    uv=[r for r in rows if r['Pollutant'].strip()=='UV-B']
    if (len(alan),len(uv))!=(114,118):raise ValueError('Source exposure counts changed')
    # Reproduction is selected by both author tier and endpoint, not effect sign.
    reproduction=[r for r in alan if r['Component'].strip()=='Pollination success']
    if any(r['Measurement'].strip() not in ENDPOINTS for r in reproduction):
        raise ValueError('Unexpected reproductive outcome requires explicit definition')
    source_effects=[{k:r[k] for k in ('NR','Study','DOI','Measurement','Species','Light.type',
        'Field.lab','SMD','SE','Mean.Treatment','SD.Treatment','N.Treatment',
        'Mean.Control','SD.Control','N.Control')} for r in reproduction]
    output={'schema':'ALAN_original_Guenat486_endpoint_paper_dependence_sensitivity_20261114_v1',
        'source_receipt':receipt,'all_source_effects':len(rows),
        'all_source_papers':len({r['Study'] for r in rows}),
        'ALAN_effects':len(alan),'ALAN_papers':len({r['Study'] for r in alan}),
        'UV_B_effects':len(uv),'UV_B_papers':len({r['Study'] for r in uv}),
        'ALAN_measurement_counts':dict(Counter(r['Measurement'] for r in alan)),
        'ALAN_reproduction':summarize(reproduction),
        'ALAN_reproduction_field_only':summarize([r for r in reproduction if r['Field.lab']=='Field']),
        'original_published_reproductive_effect_rows':source_effects,
        'limitations':[
            'Secondary extraction from published 2023 meta dataset, not original experimental plant records',
            'Uses rounded archived SMD and SE; does not reproduce original rma.mv model that recomputed some effect sizes',
            'Within-paper shared controls and treatment dependence modeled only by explicit common-correlation sensitivity',
            'Averaging distinct treatments within paper changes the estimand to a study-mean response',
            'Author seed-set/fruit-set category labels require verification of operational denominators before biological interpretation',
            'Publication identity is not proof of independently randomized biological experiments',
            'No direct ALAN-versus-UV moderator coefficient, no mixed-endpoint overall mean, no island interaction fitted',
            '2023 source inventory is not exhaustive through November 2026; no global null claim from non-significance']}
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'alan_reproductive_meta_study_level.json').write_text(json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False))
    print('ACTUAL_ALAN_REPRODUCTIVE_META',json.dumps(output,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    main(p.parse_args().out)
