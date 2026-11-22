#!/usr/bin/env python3
"""Authentic original 16-flume 2022 ALAN × crayfish emergence audit.

Source: Ogbeide et al. 2026 DOI 10.1111/1365-2435.70335 / Zenodo 10.5281/zenodo.14516568,
source repository colinog/ALAN_crayfish_SIA pinned to public commit. Experiment is
shared by Manfrin et al. 2025 DOI 10.1016/j.gecco.2025.e03577 and Arias et al. 2026
DOI 10.1186/s12302-026-01345-4. All are one 16-flume physical programme.

This script does NOT fit authors' GLMM or reconstruct spider abundance.
It downloads only small public RDS sources, verifies source Git blob SHA, uses
base-R to convert RDS into ephemeral CSV, and emits aggregated effects with
n=4 real flumes per arm. No pseudo-independent week/taxon/paper observations.
"""
from __future__ import annotations
import argparse, collections, csv, hashlib, itertools, json, math, os, pathlib, statistics, subprocess, tempfile, urllib.parse, urllib.request

SOURCE_REPO = "colinog/ALAN_crayfish_SIA"
SOURCE_COMMIT = "5e749c03161ad9ae7d57f12ddb86d9aaad283077"
SOURCES = {
 "emergence": {
  "path": "data/other data_output/emergence.rds",
  "blob_sha": "53dfb2f8940f1a6fb54ef5daf4e39699e16fb6ad",
 },
 "other_emergence": {
  "path": "data/other data_output/other_emergence.rds",
  "blob_sha": "472ec0fb0334f5db042a3a774b57ce5181fdeb36",
 },
}
TREATMENTS = ("control", "alan", "cray", "alan+cray")


def source_url(path):
 if path not in [d["path"] for d in SOURCES.values()]:
  raise ValueError("not a frozen source path")
 return f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{urllib.parse.quote(path, safe='/')}"


def git_blob_sha(data):
 return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\x00" + data).hexdigest()


def read_and_verify_rds(path, sha, limit=100000):
 url=source_url(path)
 req=urllib.request.Request(url,headers={"User-Agent":"ALAN-source-specific-science-audit/2026-11-23"})
 with urllib.request.urlopen(req,timeout=35) as resp:
  if resp.geturl().split('/')[2] != 'raw.githubusercontent.com': raise ValueError("source redirect changed hosts")
  content=resp.read(limit+1)
 if len(content)>limit:raise ValueError("source file exceeds bounded read")
 got=git_blob_sha(content)
 if got != sha:raise ValueError(f"original RDS blob SHA mismatch: {got} != {sha}")
 return content, {"original_repository":SOURCE_REPO,"original_commit":SOURCE_COMMIT,"path":path,"bytes":len(content),"git_blob_sha1":got,"sha256":hashlib.sha256(content).hexdigest(),"source_url":url}


def as_records(path):
 with open(path,newline="",encoding="utf-8-sig") as f:
  return list(csv.DictReader(f))


def source_grain(rows, *, primary=True):
    """Original Chironomidae RDS: 3 traps x 2 days per real flume-week.

    We aggregate equal-area, equal-duration sub-samples to one mean emergence
    flux per flume-week. Independent randomised lamps/flumes remain n=16.
    Do not aggregate other emergence taxa with missing taxon-specific rows.
    """
    if not rows: raise ValueError('original RDS has no records')
    expected={'treat','week','flume','abund_rate'}
    if not expected<=set(rows[0]):raise ValueError('original expected column mismatch '+str(sorted(rows[0])))
    if not primary:
        raise ValueError('Other-emergence taxa lack a balanced per-taxon collection frame; SOURCE_STOP for effect fitting')
    needed={'treat','week','flume','trap','day','Family','abundance','trap_day','area','abund_rate'}
    if not needed<=set(rows[0]):raise ValueError('original six-subsample schema changed; stop')
    values=collections.defaultdict(dict);unique_flumes={};weeks=collections.Counter();families=set()
    for r in rows:
        trt=r['treat'].strip().lower();week=r['week'].strip();flume=str(r['flume']).strip()
        if trt not in TREATMENTS:raise ValueError('Unexpected original treatment '+repr(trt))
        if week not in ('Week 1','Week 4'):raise ValueError('Unexpected original sample week '+repr(week))
        if not flume:raise ValueError('missing flume')
        if flume in unique_flumes and unique_flumes[flume]!=trt:raise ValueError('flume changes lighting assignment')
        unique_flumes[flume]=trt
        trap=r['trap'].strip();day=r['day'].strip();family=r['Family'].strip()
        families.add(family)
        if (trap,day) not in set(itertools.product(('up','mid','down'),('d1','d2'))):
            raise ValueError('original trap-day labels outside frozen six-subsample design')
        duration=float(r['trap_day']);area=float(r['area']);abundance=float(r['abundance']);rate=float(r['abund_rate'])
        if not all(map(math.isfinite,[duration,area,abundance,rate])) or duration<=0 or area<=0 or abundance<0 or rate<0:
            raise ValueError('invalid original emergence effort or count')
        if (duration,area)!=(4.0,1.0):raise ValueError('original effort or area unequal; require fresh weighted estimator')
        if not math.isclose(abundance/(duration*area),rate,rel_tol=1e-8,abs_tol=1e-8):
            raise ValueError('original rate does not equal count/effort')
        key=(trt,week,flume)
        if (trap,day) in values[key]:raise ValueError('duplicate original independent flume-week/stratum trap-day')
        values[key][(trap,day)]=rate
        weeks[week]+=1
    if families != {'chiro'}:raise ValueError('unexpected family not Chironomidae')
    if len(unique_flumes)!=16 or set(unique_flumes.values())!=set(TREATMENTS):
        raise ValueError('16 original flumes and all four treatments not confirmed')
    arm_flumes={a:{f for f,v in unique_flumes.items() if v==a} for a in TREATMENTS}
    if [len(arm_flumes[t]) for t in TREATMENTS]!=[4,4,4,4]:raise ValueError('4 flumes per arm not confirmed')
    expected_grid=set(itertools.product(('up','mid','down'),('d1','d2')))
    if len(values)!=32 or set(weeks)!={'Week 1','Week 4'}:
        raise ValueError('original 16 flumes x 2 repeated weeks incomplete')
    for arm in TREATMENTS:
        for week in ('Week 1','Week 4'):
            for flume in arm_flumes[arm]:
                key=(arm,week,flume)
                if set(values.get(key,{}))!=expected_grid:raise ValueError('unequal number of source traps/days; no pseudo replication')
    bykey={k:statistics.mean(v.values()) for k,v in values.items()}
    grain={'n_rows':len(rows),'n_real_flumes':len(unique_flumes),
           'physical_unit':'flume','n_flume_week_means':len(bykey),
           'subsamples_per_flume_week':len(expected_grid),
           'subsample_design':'3 trap positions x 2 equally-weighted days; 4 day effort, 1 area per original row',
           'treatments':dict(collections.Counter(unique_flumes.values())),
           'weeks':dict(weeks),'taxon_labels':sorted(families),
           'source_subsamples_not_independent_lamps':True}
    return bykey,grain

def log_rr(x,y):
 """ln ratio of arithmetic arm means; non-zero means only, independent physical flumes."""
 if not x or not y or min(statistics.mean(x),statistics.mean(y))<=0:raise ValueError("zero arm mean: no pseudo counts")
 mx,my=statistics.mean(x),statistics.mean(y)
 sx,sy=statistics.variance(x),statistics.variance(y)
 ln=math.log(mx/my)
 v=sx/(len(x)*mx*mx)+sy/(len(y)*my*my)
 return ln,math.sqrt(v),v


def permutation_p(x,y):
 """Exact conditional 4-of-8 permutation, no outcome-selection or pseudoreplication."""
 obs=abs(log_rr(x,y)[0]);allvals=x+y; m=len(x);n=len(allvals)
 more=0;valid=0
 for left in itertools.combinations(range(n),m):
  left=set(left); a=[z for j,z in enumerate(allvals) if j in left];b=[z for j,z in enumerate(allvals) if j not in left]
  if min(statistics.mean(a),statistics.mean(b))<=0:continue
  val=abs(log_rr(a,b)[0])
  valid+=1
  if val+1e-12>=obs:more+=1
 return more/valid if valid else None,valid


def primary_results(rows):
 bykey,grain=source_grain(rows)
 units=sorted({f for (_,_,f) in bykey})
 assignment={f:t for t,w,f in bykey}
 periods=sorted({w for t,w,f in bykey})
 estimates=[];full=[]
 for week in periods:
  cells={a:[bykey[(a,week,f)] for f in units if assignment[f]==a] for a in TREATMENTS}
  stats={a:{"n_flumes":len(v),"arithmetic_mean":statistics.mean(v),"sample_sd":statistics.stdev(v),"min":min(v),"max":max(v)}for a,v in cells.items()}
  # Explicit pre-defined factorial contrasts. Distinct physical flumes in all groups.
  contrasts=[("ALAN_no_cray","alan","control"), ("ALAN_with_cray","alan+cray","cray"),("Cray_without_ALAN","cray","control"),("Cray_with_ALAN","alan+cray","alan")]
  for contrast,a,b in contrasts:
   ln,se,v=log_rr(cells[a],cells[b]);p,k=permutation_p(cells[a],cells[b])
   estimates.append({"week":week,"contrast":contrast,"treatment":a,"reference":b,"n_treatment_flumes":len(cells[a]),"n_reference_flumes":len(cells[b]),"lnRR":ln,"RR":math.exp(ln),"delta_SE":se,"normal95_descriptive":[math.exp(ln-1.96*se),math.exp(ln+1.96*se)],"exact_4of8_label_permutation_two_sided_p":p,"permutation_allocations":k,"adopted_into_main_meta":False})
  a=log_rr(cells['alan+cray'],cells['cray'])
  b=log_rr(cells['alan'],cells['control'])
  dint=a[0]-b[0];se=math.sqrt(a[2]+b[2])
  estimates.append({"week":week,"contrast":"factorial_ALAN_x_Cray_interaction","n_independent_flumes":16,"lnRR":dint,"RR":math.exp(dint),"delta_SE":se,"normal95_descriptive":[math.exp(dint-1.96*se),math.exp(dint+1.96*se)],"exact_permutation_p":None,"adopted_into_main_meta":False})
  full.append({"week":week,"summary":stats})
 # Within-flume time comparison: individual log flux if zero? explicitly separate; do not inflate k.
 return {"grain":grain,"unit_counts":full,"contrasts":estimates,"warning":"Two weeks from same 16 flumes are dependent, 4 arms are same factorial programme; n=4 per arm. Ratios/normal approximations descriptive; original author model used Gamma GLMM with AR(1)."}


def primary_fullrun(output):
 import shutil
 if not shutil.which('Rscript'):raise RuntimeError('Rscript unavailable; do not guess RDS values')
 receipts=[]
 with tempfile.TemporaryDirectory(prefix='alan-landau16-') as td:
  rows=[];other=[]
  for k,d in SOURCES.items():
   content,proof=read_and_verify_rds(d['path'],d['blob_sha'])
   receipts.append(proof)
   rds=pathlib.Path(td)/(k+'.rds');csvpath=pathlib.Path(td)/(k+'.csv')
   rds.write_bytes(content)
   script=pathlib.Path(__file__).with_name('alan_landau16_rds_to_csv_20261123.R')
   cp=subprocess.run(['Rscript',str(script),str(rds),str(csvpath)],capture_output=True,text=True,timeout=60)
   if cp.returncode:raise RuntimeError('RDS original reader failed: '+cp.stderr[-1200:])
   if k=='emergence':rows=as_records(csvpath)
   else:other=as_records(csvpath)
  res=primary_results(rows)
  # Secondary source contains emergence families, but not equivalent to main and never pooled.
  og=original_grain_diagnostic(other)
  og['source_grain_status']='HETEROGENEOUS_TAXA_NOT_COMBINED'
  og['unbalanced_taxon_rows_preclude_straight_sum']=True
  result={"study":"Ogbeide et al. 2026 / Manfrin et al. 2025 shared physical 2022 Landau RSM programme", "programme_id":"Landau_RSM_2022_16flumes_ALAN_crayfish","doi":["10.1016/j.gecco.2025.e03577","10.1111/1365-2435.70335","10.1186/s12302-026-01345-4"],"authentic_source_doi":"10.5281/zenodo.14516568","original_source_receipts":receipts,"source_reference_code":"https://github.com/colinog/ALAN_crayfish_SIA/blob/5e749c03161ad9ae7d57f12ddb86d9aaad283077/code/other_code_output.R", "emergence_dataset":"Chironomidae emergence flux (author emergence.rds, not air-eclector capture)","other_emergence_grain":og,"original_primary":res,"number_new_independent_programmes":0,"source_independence_note":"Already entered under Ogbeide2026 in ALAN studies.csv; distinct response study, no addition to existing 108 effect catalogue until adjudication", "inference":"Exploratory source-only independent 4-flume-per-arm arithmetic mean lnRR; do not interpret as main wet-plant mechanism or causal spider-mediated effect."}
  pathlib.Path(output).parent.mkdir(parents=True,exist_ok=True)
  pathlib.Path(output).write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
  print('ORIGINAL_LANDAU_FLUME_RESULTS',len(res['contrasts']),'weeks',res['grain']['weeks'],'RESULT',output)
  for x in res['contrasts']:
   print('CONTRAST',x['week'],x['contrast'],'RR',round(x['RR'],4),'n',x.get('n_treatment_flumes',4),'p',x.get('exact_4of8_label_permutation_two_sided_p'))
  return result


def original_grain_diagnostic(rows):
    if not rows: raise ValueError('Original RDS empty')
    columns=list(rows[0].keys())
    summary={}
    for col in columns:
        nonempty=[str(r.get(col, '')).strip() for r in rows if str(r.get(col,'')).strip()]
        uniq=sorted(set(nonempty))
        summary[col]={'nonempty':len(nonempty),'distinct':len(uniq)}
        if len(uniq)<=25:summary[col]['values']=uniq
    group=collections.Counter((r.get('treat','').strip(),r.get('week','').strip(),r.get('flume','').strip()) for r in rows)
    flume_treat=collections.defaultdict(set)
    for r in rows:flume_treat[r.get('flume','').strip()].add(r.get('treat','').strip())
    return {'n_rows':len(rows),'columns':summary,'n_flume_week_groups':len(group),
        'rows_per_flume_week_histogram':dict(collections.Counter(group.values())),
        'treatment_flume_pairs':len({(r.get('treat','').strip(),r.get('flume','').strip()) for r in rows}),
        'flumes_assigned_multiple_treatments':{k:sorted(v) for k,v in flume_treat.items() if len(v)>1},
        'sample_group_counts':{str(k):v for k,v in sorted(group.items())[:20]}}


def inspect_original_grain(output):
    import shutil
    if not shutil.which('Rscript'):raise RuntimeError('Rscript unavailable')
    receipts=[];data={}
    with tempfile.TemporaryDirectory(prefix='alan-landau-inspection-') as td:
      for key,d in SOURCES.items():
        raw,proof=read_and_verify_rds(d['path'],d['blob_sha']);receipts.append(proof)
        rds=pathlib.Path(td)/(key+'.rds');csvpath=pathlib.Path(td)/(key+'.csv');rds.write_bytes(raw)
        cp=subprocess.run(['Rscript',str(pathlib.Path(__file__).with_name('alan_landau16_rds_to_csv_20261123.R')),str(rds),str(csvpath)],capture_output=True,text=True,timeout=60)
        if cp.returncode:raise RuntimeError('RDS read error '+cp.stderr[-500:])
        data[key]=original_grain_diagnostic(as_records(csvpath))
    result={'status':'SOURCE_GRAIN_AUDIT_ONLY_NO_EFFECT','sources':receipts,'dataset_grain':data,
         'no_effects_estimated':True,'no_new_programme_k':True}
    pathlib.Path(output).parent.mkdir(parents=True,exist_ok=True)
    pathlib.Path(output).write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('GRAIN_AUDIT',json.dumps({k:{'n_rows':v['n_rows'],'rows_per_flume_week_histogram':v['rows_per_flume_week_histogram'], 'columns':list(v['columns'])} for k,v in data.items()}))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',default='results/alan_landau16_original_rds_20261123.json')
    p.add_argument('--inspect-grain',action='store_true')
    args=p.parse_args()
    if args.inspect_grain:inspect_original_grain(args.output)
    else:primary_fullrun(args.output)

if __name__=='__main__':main()
