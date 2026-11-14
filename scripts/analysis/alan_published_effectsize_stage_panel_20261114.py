"""Stage-specific ALAN summary-data reanalysis. No causal pollination claim.
Original byte-identical Sanders archive; raw ratio signs, not ecosystem value.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import t


def log_rr(mt, sdt, nt, mc, sdc, nc):
    a=np.asarray([mt,sdt,nt,mc,sdc,nc],dtype=float)
    if not np.isfinite(a).all() or mt<=0 or mc<=0 or sdt<0 or sdc<0 or nt<2 or nc<2:
        raise ValueError('Finite positive means, nonnegative SD, n >= 2 required')
    v=sdt*sdt/(nt*mt*mt)+sdc*sdc/(nc*mc*mc)
    if v<=0:raise ValueError('Positive sampling variance required')
    return math.log(mt/mc),v


def bounded_sd_ok(mean,sd,n,lower,upper):
    if not all(math.isfinite(x) for x in (mean,sd,n,lower,upper)) or n<2 or sd<0 or upper<=lower:
        return False
    if not lower<=mean<=upper:return False
    return sd*sd <= n/(n-1)*(mean-lower)*(upper-mean)+1e-10*max(1,upper*upper)


def nonnegative_sd_ok(mean,sd,n):
    # For ordinary nonnegative observations, sample SD <= sqrt(n)*sample mean.
    if not all(math.isfinite(x) for x in (mean,sd,n)) or mean<0 or sd<0 or n<2:return False
    return sd <= math.sqrt(n)*mean+1e-10*max(1.,mean)


def aggregate(ys,vs,rho):
    y=np.asarray(ys,dtype=float);v=np.asarray(vs,dtype=float)
    if len(y)==0 or len(y)!=len(v) or not 0<=rho<=1 or not np.isfinite(y).all() or not np.isfinite(v).all() or np.any(v<=0):
        raise ValueError('Matching finite effects/positive variances and rho in [0,1] required')
    m=len(y)
    var=(float(v.sum())+rho*(float(np.sqrt(v).sum())**2-float(v.sum())))/(m*m)
    return float(y.mean()),var


def re_meta(ys,vs):
    y=np.asarray(ys,dtype=float);v=np.asarray(vs,dtype=float);k=len(y)
    if k<2 or len(v)!=k or not np.isfinite(y).all() or not np.isfinite(v).all() or np.any(v<=0):
        raise ValueError('At least two finite study estimates and positive variances required')
    def objective(tau2):
        w=1/(v+tau2);mu=float(np.sum(w*y)/np.sum(w))
        return .5*(np.log(v+tau2).sum()+math.log(float(w.sum()))+np.sum(w*(y-mu)**2))
    upper=max(1.,float(np.var(y,ddof=1))*100.,float(v.max())*10.)
    fit=minimize_scalar(objective,bounds=(0.,upper),method='bounded',options={'xatol':1e-12})
    if not fit.success:raise ValueError('REML optimization failed')
    tau2=0. if objective(0)<=fit.fun else float(fit.x)
    w=1/(v+tau2);mu=float(np.sum(w*y)/np.sum(w));q=float(np.sum(w*(y-mu)**2))
    scale=max(1.,q/(k-1));se=math.sqrt(scale/float(w.sum()));crit=float(t.ppf(.975,k-1))
    ci=[mu-crit*se,mu+crit*se]
    return {'k_studies':k,'estimate':mu,'tau2':tau2,'se_mkh':se,'q':q,'kh_scale':scale,'df':k-1,
            'ci95':ci,'ratio':math.exp(mu),'ratio_ci95':[math.exp(a) for a in ci],
            'p_two_sided':float(2*t.sf(abs(mu/se),k-1)),
            'method':'REML; modified Hartung-Knapp scale max(1,Q/(k-1)); t(k-1)'}


CAPTURE_ENDPOINTS={
 'Handled.prey.(No.predatitor.cue)','Handled.prey.(With.predator.cue)',
 'Consumption.of.insects','prey/min','fish.consumption.sandy.base',
 'fish.consumption','small.prey.captured'}
LABELS={'2':'Underwood 2017 / intertidal gastropod',
        '160':'Rapatsa & Moyo 2017 / aquaculture fish',
        '387':'Santos 2010 / shorebirds',
        '501':'Czarnecka 2019 / perch',
        '525':'Willmott 2019 / orb-web spiders'}
SOURCE_MD5='665c1077ea63b2e11ededd24371ce0cf'
NUMERIC=('Experimental_Mean','Experimental_SD','Experimental_N','Control_Mean','Control_SD','Control_N')


def read_source(path):
    import csv,hashlib,io
    from pathlib import Path
    b=Path(path).read_bytes()
    if hashlib.md5(b).hexdigest()!=SOURCE_MD5:raise ValueError('Source snapshot MD5 mismatch')
    raw_rows=list(csv.DictReader(io.StringIO(b.decode('latin1'),newline=None)))
    rows=[{k.strip():v.strip() for k,v in r.items()} for r in raw_rows]
    if len(rows)!=1304 or len(rows[0])!=30:raise ValueError('Original source dimensions changed')
    return rows,{'archive':'Zenodo 3998338; 2020 snapshot, not verified latest Dryad revision',
                 'dataset_doi':'10.5061/dryad.wpzgmsbjn','article_doi':'10.1038/s41559-020-01322-x',
                 'md5':SOURCE_MD5,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b),
                 'decoding':'latin1 reversible, universal newlines; original bytes untouched'}


def source_contrast(r):
    a={k:float(r[k]) for k in NUMERIC}
    y,v=log_rr(*[a[k] for k in NUMERIC]);flags=[]
    for arm in ('Control','Experimental'):
        if not nonnegative_sd_ok(a[arm+'_Mean'],a[arm+'_SD'],a[arm+'_N']):
            flags.append(arm+'_SD_incompatible_with_nonnegative_sample_mean_and_N')
    if r['Study']=='160':
        flags.append('source_N100_not_verified_independent_light_units_primary_paper_three_tanks_per_treatment')
    return {'study':r['Study'],'source_id':r['ID'],'endpoint':r['Measure.description'],
            'doi':r['DOI'],'species':r['Species'],'log_response_ratio':y,'variance_from_source_summaries':v,
            'ratio':math.exp(y),'source_numeric':a,'flags':flags}


def synthesize(contrasts,rho):
    from collections import defaultdict
    grouped=defaultdict(list)
    for r in contrasts:grouped[r['study']].append(r)
    studies=[]
    for study in sorted(grouped,key=int):
        rr=grouped[study]
        y,v=aggregate([a['log_response_ratio'] for a in rr],[a['variance_from_source_summaries'] for a in rr],rho)
        studies.append({'study':study,'label':LABELS[study],'n_contrasts':len(rr),
                        'log_response_ratio':y,'within_study_variance':v,'geometric_mean_ratio':math.exp(y)})
    result=re_meta([r['log_response_ratio'] for r in studies],[r['within_study_variance'] for r in studies])
    result.update({'rho_within_study':rho,'n_contrasts':len(contrasts),'studies':studies,
                   'causal_pooled_moth_predation_effect':False,
                   'study_id_is_publication_cluster_not_a_verified_independent_light_installation':True})
    w=np.array([1/(r['within_study_variance']+result['tau2']) for r in studies])
    for r,weight in zip(studies,w/w.sum()):r['normalized_meta_weight']=float(weight)
    if len(studies)>2:
        result['leave_one_publication_out']=[]
        for omit in studies:
            remain=[r for r in studies if r['study']!=omit['study']]
            z=re_meta([r['log_response_ratio'] for r in remain],[r['within_study_variance'] for r in remain])
            result['leave_one_publication_out'].append({'omitted_study':omit['study'],
                'ratio':z['ratio'],'ratio_ci95':z['ratio_ci95'],'p_two_sided':z['p_two_sided']})
    return result


def run(source,outdir):
    import collections,csv,json
    from pathlib import Path
    rows,receipt=read_source(source);out=Path(outdir)
    if out.exists() and any(out.iterdir()):raise ValueError('Refuse overwrite of result directory')
    out.mkdir(parents=True,exist_ok=True)
    kept=[r for r in rows if r['Delete']=='Keep']
    selected=[source_contrast(r) for r in kept if r['Measure.description'] in CAPTURE_ENDPOINTS]
    if len(selected)!=14 or {r['study'] for r in selected}!=set(LABELS):raise ValueError('Source support changed')
    math_screened=[r for r in selected if not any('SD_incompatible' in f for f in r['flags'])]
    guarded=[r for r in selected if not r['flags']]
    scenarios={name:[synthesize(rr,rho) for rho in (0.,.5,1.)] for name,rr in
               [('unfiltered_import_diagnostic',selected),('nonnegative_SD_screen_only',math_screened),
                ('SD_and_known_replication_guard',guarded)]}
    plant=[]
    for r in kept:
        if r['Study']=='68' and r['Measure.description']=='pollen.presence':
            z=source_contrast(r);z['stage']='moth_carried_pollen';z['primary_verified_doi']='10.1111/gcb.13371'
            z['source_doi_year_drift_warning']=True;plant.append(z)
        if r['Study']=='72' and r['Measure.description'] in {'Number.of.Visits','Number.of.fruits','Percentage.of.fruits'}:
            z=source_contrast(r);z['stage']='visits' if r['Measure.description']=='Number.of.Visits' else 'plant_reproduction';plant.append(z)
        if r['Study']=='538' and r['Measure.description']=='pollination.success':
            z=source_contrast(r);z['stage']='plant_pollination_success';z['lighting_arm_mapping_unresolved']=True
            z['bounded_percentage_SD_checks']={arm:bounded_sd_ok(float(r[arm+'_Mean']),float(r[arm+'_SD']),float(r[arm+'_N']),0.,100.) for arm in ('Control','Experimental')}
            z['exclude_from_inverse_variance_plant_synthesis']=not all(z['bounded_percentage_SD_checks'].values())
            plant.append(z)
    if len([p for p in plant if p['study']=='538'])!=4:raise ValueError('Source pollination subset changed')
    repro=[r for r in kept if r['Sub.Category']=='Reproductive.output']
    flowering=[r for r in kept if r['Measure.description'].startswith('Flowering.date')]
    bystudy=collections.defaultdict(list)
    for r in rows:bystudy[r['Study']].append(r)
    conflicts=[]
    for sid,rr in bystudy.items():
        ds=sorted(set(r['DOI'] for r in rr));ys=sorted(set(r['Publication.year'] for r in rr))
        if len(ds)>1 or len(ys)>1:conflicts.append({'study':sid,'distinct_doi':len(ds),'distinct_years':ys})
    counts={'source_rows':len(rows),'source_study_ids':len(bystudy),'source_doi_strings':len(set(r['DOI'] for r in rows)),
            'author_Keep_rows':len(kept),'author_Delete_rows':len(rows)-len(kept),
            'rows_year_after_archive_2020':sum(float(r['Publication.year'])>2020 for r in rows),
            'reproductive_output_rows':len(repro),'reproductive_output_studies':len(set(r['Study'] for r in repro)),
            'known_plant_reproduction_study_ids':['72','538'],
            'flowering_date_rows':len(flowering),'flowering_date_study_ids':sorted(set(r['Study'] for r in flowering))}
    result={'schema':'ALAN_published_summary_stage_panel_v1','source':receipt,'counts':counts,
            'metadata_conflicts':conflicts,'predation_selected':selected,'predation_scenarios':scenarios,
            'plant_function_contrasts':plant,
            'not_fitted':{'global_plant_meta_effect':True,'phenology_mismatch':True,
                          'island_mainland_interaction':True,'predator_pollination_compensation':True},
            'interpretation_limits':[
             'New calculations from published extracted summaries, not raw primary observations or reproduction of original SMDH/MCMC model.',
             'Raw log response ratio signs indicate more/less prey handling/capture/consumption, not ecological benefit or damage.',
             'Capture panel mixes intertidal, freshwater and terrestrial consumers; not adult-moth-predation or plant-pollination.',
             'Known numeric/replication issues screened; remaining SD/N and paired designs not exhaustively primary-table audited.',
             'Within-paper covariance unknown; rho=0/.5/1 sensitivities, not estimated correlations.',
             'Only three publication clusters remain after both guards; no high-dimensional moderator or causal inference.',
             'Study538 impossible literal percentage SD does not invalidate its original paper; no corrected SD imputed.',
             'Source years and DOI strings have internal drift; neither used as independent-study identifiers or temporal trends.',
             'No direct predator-as-alternative-pollinator comparison observed.']}
    out.joinpath('stage_meta_result.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    with out.joinpath('capture_summary_scenarios.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['scenario','rho','n_publications','n_contrasts','ratio','ci_low','ci_high','p','tau2']);writer.writeheader()
        for name,ss in scenarios.items():
            for z in ss:writer.writerow({'scenario':name,'rho':z['rho_within_study'],'n_publications':z['k_studies'],'n_contrasts':z['n_contrasts'],'ratio':z['ratio'],'ci_low':z['ratio_ci95'][0],'ci_high':z['ratio_ci95'][1],'p':z['p_two_sided'],'tau2':z['tau2']})
    print('STAGE_SOURCE_COUNTS',json.dumps(counts))
    print('STAGE_PRIMARY_GUARDED_META',json.dumps(scenarios['SD_and_known_replication_guard'][1]))
    print('STAGE_INPUT_QUALITY_GATES',json.dumps({'quarantined_capture_rows':[(r['source_id'],r['flags']) for r in selected if r['flags']],
          'study538_percent_rows_excluded':sum(p.get('exclude_from_inverse_variance_plant_synthesis',False) for p in plant)}))
    return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.source,a.out)
