"""Replay original uploaded Dryad workbook; standard scientific analysis of numeric tables.
python scripts/analysis/alan_eeres_arrangement_verified_20261125.py DISCONNECTED_alldata.xlsx --out /tmp/arrangement.json
Install openpyxl pandas numpy scipy for local repository replay.
"""
import argparse, hashlib, json, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t
SHA="496dcfa07436619084dc677317525dc674775e9965a0fdfb482cfa1cb4d8977f"
T=t.ppf(.975,5)
def contrast(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    assert len(a)==len(b)==6
    d=a-b; m=d.mean(); s=d.std(ddof=1)/np.sqrt(6)
    assert a.mean()>0 and b.mean()>0
    th=np.log(a.mean()/b.mean());se=(a/a.mean()-b/b.mean()).std(ddof=1)/np.sqrt(6)
    return dict(num_mean=float(a.mean()),den_mean=float(b.mean()),diff=float(m),diff_ci=[float(m-T*s),float(m+T*s)],block_diffs=d.tolist(),rr=float(np.exp(th)),ci=np.exp(th+np.array([-1,1])*T*se).tolist(),se=float(se),p=float(2*t.sf(abs(th/se),5)),signflip_reference_p=sum(abs(np.mean(d*np.array(s)))>=abs(m)-1e-10 for s in itertools.product([-1,1],repeat=6))/64,leave_one_out=[float(np.delete(a,i).mean()/np.delete(b,i).mean()) for i in range(6)])
def absolute(a):
    a=np.asarray(a);m=a.mean();s=a.std(ddof=1)/np.sqrt(6)
    return dict(mean=float(m),ci=[float(m-T*s),float(m+T*s)],blocks=a.tolist())
def analyze(d,b):
    assert len(d)==144 and len(b)==108
    for frame in (d,b):
        assert frame.ID.is_unique and set(frame.HT)=={f'ht{i}' for i in range(1,7)} and frame.flume.nunique()==18
    cols=list(d.columns[8:]);bcols=list(b.columns[6:])
    assert len(cols)==38 and len(bcols)==40
    for frame,cs,total in ((d,cols,5564),(b,bcols,11326)):
        v=frame[cs].to_numpy(float)
        assert np.isfinite(v).all() and (v>=0).all() and (v%1==0).all() and v.sum()==total
        frame['count']=v.sum(1);frame['rich']=(v>0).sum(1)
    keys=['HT','flume','time'];mapping=d[['HT','flume','treatment_down']].drop_duplicates()
    assert len(mapping)==18 and (mapping.groupby(['HT','treatment_down']).size()==1).all()
    expected={'CT':['CT','CT'],'SS':['F','L'],'CS':['CT','L+F']}
    assert set(d.time)=={'T1','T2','T3','T4'}
    for _,z in d.groupby(keys):
        assert len(z)==2 and z.side.nunique()==2 and z.treatment_down.nunique()==1 and sorted(z.treatment_up)==expected[z.treatment_down.iloc[0]]
    assert len(d.groupby(keys))==72
    agg=d.groupby(keys+['treatment_down'])[cols].sum()
    events=agg.assign(count=agg.sum(axis=1),rich=(agg>0).sum(axis=1)).reset_index()
    assert events['count'].sum()==5564 and (events.groupby(['HT','flume','treatment_down']).size()==4).all()
    drift=events.groupby(['HT','treatment_down'])[['count','rich']].mean().reset_index().rename(columns={'treatment_down':'arm'})
    bd=b[b.section=='down'].copy();bd['arm']=bd.treatment.replace({'CT_down':'CT'})
    assert len(bd)==36 and (bd.groupby(['HT','flume','arm']).size()==2).all()
    check=bd.merge(mapping,on=['HT','flume'],validate='many_to_one');assert (check.arm==check.treatment_down).all()
    benthos=bd.groupby(['HT','arm'])[['count','rich']].mean().reset_index()
    out={'source_sha256':SHA,'blocks':{},'contrasts':{},'by_week':{},'decomposition':{},'stage_discordance':{}}
    for stage,frame in [('drift',drift),('down_benthos',benthos)]:
        for metric in ['count','rich']:
            p=frame.pivot(index='HT',columns='arm',values=metric)
            assert p.shape==(6,3) and p.notna().all().all()
            k=stage+'_'+metric;out['blocks'][k]={a:p[a].tolist() for a in ['CT','SS','CS']}
            out['contrasts'][k]={f'{a}/{c}':contrast(p[a],p[c]) for a,c in [('CS','SS'),('SS','CT'),('CS','CT')]}
    q=0
    for i,(p,k) in enumerate(sorted((v['CS/SS']['p'],k) for k,v in out['contrasts'].items())):
        q=max(q,min(1,p*(4-i)));out['contrasts'][k]['CS/SS']['holm4']=q
    for ti,g in events.groupby('time'):
        p=g.pivot(index='HT',columns='treatment_down',values='count')
        r=contrast(p.CS,p.SS);out['by_week'][ti]={'means':p.mean().to_dict(),'CS_SS':r['rr'],'ci':r['ci']}
    base=d[(d.treatment_down=='CT')&(d.treatment_up=='CT')].groupby('HT')['count'].mean()
    x=d[(d.treatment_down=='CS')&(d.treatment_up=='CT')].groupby('HT')['count'].mean()
    l=d[d.treatment_up=='L'].groupby('HT')['count'].mean()
    f=d[d.treatment_up=='F'].groupby('HT')['count'].mean()
    lf=d[d.treatment_up=='L+F'].groupby('HT')['count'].mean()
    p=drift.pivot(index='HT',columns='arm',values='count')
    assert np.allclose(p.CT,2*base) and np.allclose(p.CS,x+lf) and np.allclose(p.SS,l+f)
    loc=(lf-f)-(l-base);ctx=x-base;assert np.allclose(loc+ctx,p.CS-p.SS)
    br={'CT_network_per_branch':base,'CS_untreated':x,'L_SS':l,'F_SS':f,'LF_CS':lf}
    out['decomposition']=dict(branch_blocks={k:v.tolist() for k,v in br.items()},branch_means={k:float(v.mean()) for k,v in br.items()},actual_CS_untreated_fraction_of_all_CT=float(x.sum()/(2*base+x).sum()),actual_CS_untreated_vs_all_CT_mean=float(3*x.sum()/(2*base+x).sum()),untreated_branch_context_rr=contrast(x,base),local_additive_interaction_contrast=absolute(loc),untreated_branch_context_contrast=absolute(ctx),sum_equals_CS_minus_SS=absolute(loc+ctx))
    for metric in ['count','rich']:
        a=out['blocks']['drift_'+metric];b=out['blocks']['down_benthos_'+metric]
        x,y,u,v=map(np.array,[a['CS'],a['SS'],b['CS'],b['SS']])
        th=np.log(x.mean()/y.mean())-np.log(u.mean()/v.mean())
        se=(x/x.mean()-y/y.mean()-u/u.mean()+v/v.mean()).std(ddof=1)/np.sqrt(6)
        out['stage_discordance'][metric]=dict(ratio_of_RRs=float(np.exp(th)),ci=np.exp(th+np.array([-1,1])*T*se).tolist(),se=float(se))
    p=events.pivot(index='HT',columns=['time','treatment_down'],values='count')
    x,y,u,v=[p[k].to_numpy(float) for k in [('T4','CS'),('T4','SS'),('T1','CS'),('T1','SS')]]
    th=np.log(x.mean()/y.mean())-np.log(u.mean()/v.mean());se=(x/x.mean()-y/y.mean()-u/u.mean()+v/v.mean()).std(ddof=1)/np.sqrt(6)
    out['CS_SS_T4_vs_T1']=dict(rrr=float(np.exp(th)),ci=np.exp(th+np.array([-1,1])*T*se).tolist())
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('xlsx',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert hashlib.sha256(a.xlsx.read_bytes()).hexdigest()==SHA
    r=analyze(pd.read_excel(a.xlsx,sheet_name='drift'),pd.read_excel(a.xlsx,sheet_name='benthos'))
    assert abs(r['contrasts']['drift_count']['CS/SS']['rr']-0.8273743016759777)<1e-12
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2)+'\n')
