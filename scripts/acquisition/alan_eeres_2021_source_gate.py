#!/usr/bin/env python3
"""Bounded official Dryad audit of Burgazzi EERES-2021 DISCONNECTED 18-flume dataset.
No raw observations are re-published. No treatment effect is inferred at this stage.
"""
import csv, io, json, hashlib, os, re, zipfile, argparse, math, statistics
from collections import Counter
from pathlib import Path
from urllib.parse import urljoin, urlparse, quote
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

DOI='10.5061/dryad.rfj6q57k5'; B='https://datadryad.org/api/v2'
NAME='DISCONNECTED_alldata.xlsx'; MAX=2_000_000
SOURCE_SHA='496dcfa07436619084dc677317525dc674775e9965a0fdfb482cfa1cb4d8977f'
HOSTS={'datadryad.org','www.datadryad.org'}
# Exact delivery host observed after following official Dryad API version ZIP.
DELIVERY_HOST='w27i4vubmmwkssdf6xpszafmhq0rffog.lambda-url.us-west-2.on.aws'

def grab(url):
    u=urlparse(url)
    if u.scheme!='https' or u.hostname not in HOSTS:raise ValueError('UNAPPROVED_HOST')
    with urlopen(Request(url,headers={'User-Agent':'ALAN-original-source-audit/2026-11-23','Accept':'application/json,*/*'}),timeout=28) as h:
        final=h.geturl();u=urlparse(final)
        if u.scheme!='https' or u.hostname not in (HOSTS | {DELIVERY_HOST}):raise ValueError('REDIRECTED_OUTSIDE_ALLOWLIST:'+str(u.hostname))
        if int(h.headers.get('Content-Length') or '0')>MAX:raise ValueError('PAYLOAD_TOO_LARGE')
        data=h.read(MAX+1)
        if len(data)>MAX:raise ValueError('PAYLOAD_TOO_LARGE')
        return data

def hrefs(obj):
    out=set()
    def w(o):
        if isinstance(o,dict):
            for k,v in o.items():
                if k in {'href','url','link'} and isinstance(v,str):
                    u=urljoin('https://datadryad.org',v);p=urlparse(u)
                    if p.scheme=='https' and p.hostname in HOSTS and p.path.startswith('/api/v2/'):out.add(u)
                elif isinstance(v,(dict,list)):w(v)
        elif isinstance(o,list):
            for v in o:w(v)
    w(obj);return out

def inspect_xlsx(content):
    """Source grain only: calculate sheet counts, headers and categorical keys."""
    try:
        import openpyxl
    except ImportError:return {'status':'SOURCE_STOP_MISSING_OPENPYXL'}
    book=openpyxl.load_workbook(io.BytesIO(content),read_only=True,data_only=True)
    report={}
    for name in book.sheetnames:
        ws=book[name]
        itr=ws.iter_rows(values_only=True)
        try:header=next(itr)
        except StopIteration: report[name]={'status':'EMPTY_SHEET'};continue
        cols=[str(v or '').strip() for v in header]
        field_n=len(cols);freq={k:Counter() for k in cols[:min(15,len(cols))] if re.search('ht|flume|treatment|side|section|time|date|flow|light|group|block',k,re.I)}
        n=0;sum_taxa=0;negative=0;na=0
        for row in itr:
            if not any(v is not None and str(v).strip() for v in row):continue
            n+=1
            for j,k in enumerate(cols[:min(15,len(cols))]):
                if k in freq and j<len(row):
                    v=row[j]
                    if v is not None and str(v).strip()!='':freq[k][str(v)[:65]]+=1
            if n>20000:raise ValueError('TOO_MANY_ROWS')
        report[name]={'n_nonempty_data_rows':n,'n_columns':field_n,'columns':cols[:80],
                      'physical_design_candidate_columns':{k:{'distinct':len(v),'top':v.most_common(20) if len(v)<=25 else v.most_common(6)} for k,v in freq.items()}}
    return {'status':'SOURCE_EXCEL_SHEETS_READ','sheets':report}

def source_effects(content):
    """Exploratory physical-HT source comparisons; no raw record export.

    Six physical header-tank blocks, three flumes per block. One L and F
    upstream section and one L+F upstream section per block; three CT sections.
    Sample nights / benthos subsamples are *within-block* measurements.
    """
    import openpyxl
    wb=openpyxl.load_workbook(io.BytesIO(content),read_only=True,data_only=True)
    blocks={}
    for sheet_name,expected_treatment,metadata_count,counts in (
            ('drift','treatment_up',8,{'CT':12,'L':4,'F':4,'L+F':4}),
            ('benthos','treatment',6,{'CT':6,'L':2,'F':2,'L+F':2})):
        ws=wb[sheet_name]
        itr=ws.iter_rows(values_only=True)
        headers=[str(x or '') for x in next(itr)]
        if headers[1:3]!=['HT','flume'] or expected_treatment not in headers:
            return {'status':'SOURCE_STOP_SCHEMA_CHANGED','sheet':sheet_name}
        idx={k:headers.index(k) for k in ['HT','flume',expected_treatment]}
        abundance={};origin=set();invalid=0;observed_rows=0
        for row in itr:
            if not row or row[0] is None:continue
            if sheet_name=='benthos' and str(row[headers.index('section')]).strip()!='up':continue
            tr=str(row[idx[expected_treatment]]).strip()
            if tr=='CT_up':tr='CT'
            if tr not in counts:continue
            block=str(row[idx['HT']]).strip(); flume=str(row[idx['flume']]).strip()
            if block not in {f'ht{i}' for i in range(1,7)}:return {'status':'SOURCE_STOP_BLOCK_ID'}
            origin.add(flume)
            taxon=row[metadata_count:len(headers)]
            if len(taxon)!=len(headers)-metadata_count:return {'status':'SOURCE_STOP_COLUMN_MISMATCH'}
            if any(v is None or isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(float(v)) or v<0 for v in taxon):
                invalid+=1;continue
            val=sum(taxon); rich=sum(v>0 for v in taxon)
            entry=abundance.setdefault((block,tr),{'abund':[],'rich':[]})
            entry['abund'].append(val); entry['rich'].append(rich)
            observed_rows+=1
        if invalid:return {'status':'SOURCE_STOP_MISSING_OR_NONNUMERIC_TAXA','sheet':sheet_name,'n_invalid_rows':invalid}
        if len(origin)!=18:return {'status':'SOURCE_STOP_FLUME_SUPPORT','sheet':sheet_name,'flumes':len(origin)}
        per_block={};replicates={}
        for i in range(1,7):
            block=f'ht{i}'
            per_block[block]={}
            for tr,expected in counts.items():
                ent=abundance.get((block,tr))
                if not ent or len(ent['abund'])!=expected:
                    return {'status':'SOURCE_STOP_UNBALANCED_TREATMENT','sheet':sheet_name,'block':block,
                            'treatment':tr,'found':len(ent['abund']) if ent else 0,'expected':expected}
                per_block[block][tr]={k:statistics.mean(v) for k,v in ent.items()}
                replicates[f'{block}|{tr}']=len(ent['abund'])
        blocks[sheet_name]={'rows':observed_rows,'n_physical_flumes':len(origin),
                             'n_independent_header_blocks':6,'reps':replicates,'by_block':per_block}
    tcrit=2.570581835636305 # two-sided t 0.975, df=5
    outcomes={}
    for sheet_name in ('drift','benthos'):
        outcomes[sheet_name]={}
        for measure in ('abund','rich'):
            out=[]
            for tr,ref in [('L','CT'),('L+F','F')]:
                pairs=[(blocks[sheet_name]['by_block'][f'ht{i}'][tr][measure],
                        blocks[sheet_name]['by_block'][f'ht{i}'][ref][measure]) for i in range(1,7)]
                a=[x for x,y in pairs];b=[y for x,y in pairs]
                ma=statistics.mean(a); mb=statistics.mean(b)
                if ma<=0 or mb<=0:return {'status':'SOURCE_STOP_NONPOSITIVE_ARM_MEAN','sheet':sheet_name,'measure':measure}
                cov=statistics.covariance(a,b)
                va=statistics.variance(a);vb=statistics.variance(b)
                n=6
                v=(va/ma**2+vb/mb**2-2*cov/(ma*mb))/n
                if v< -1e-9:return {'status':'SOURCE_STOP_INVALID_COVARIANCE','sheet':sheet_name,'measure':measure,'v':v}
                se=math.sqrt(max(v,0.0)); logrr=math.log(ma/mb)
                out.append({'treatment':tr,'control':ref,'mean_treatment':ma,'mean_control':mb,
                            'RR':ma/mb,'lnRR':logrr,'SE_paired_HT_delta':se,
                            't5_exploratory_95RR':[math.exp(logrr-tcrit*se),math.exp(logrr+tcrit*se)],
                            'physical_HT_n':n,'block_treatment_values':a,'block_control_values':b,
                            'status':'POSTHOC_DESCRIPTIVE_NOT_META_ADOPTED'})
            outcomes[sheet_name][measure]=out
    return {'status':'SOURCE_DERIVED_HT_BLOCK_COMPARISONS','source_sha256':SOURCE_SHA,
            'physical_programme':'EERES2021_Disconnected18flumes_6headers',
            'source_grain':{'drift':'6 HT x 4 treatment classes; 4 repeated weeks; 144 rows',
                            'benthos_up':'6 HT x 4 treatment classes; 2 subsamples per section; 72 source rows'},
            'n_independent_header_blocks':6,
            'effect_summaries':outcomes,
            'joint_shared_HT':joint_source_estimands(outcomes),
            'critical_boundary':'HT grouping avoids 18-flume or repeated-week pseudoreplication. Drift is in-stream transport, not aquatic emergence to terrestrial predators; benthos a separate terminal sample. No 2026/2024 paper double count.'}


def joint_source_estimands(outcomes):
    """A shared-6-header-tank delta covariance for all 8 source log ratios.

    Distinct measurement endpoints share headers; all simultaneous 2021 light
    treatments use the same original header-based assignment. Returned results
    are post-hoc descriptors with df5 approximate intervals, not a meta sample.
    """
    if set(outcomes) != {'drift', 'benthos'}:
        return {'status':'SOURCE_STOP_MISSING_STRATA'}
    names=[]; gradients=[]; blockcolumns={}; means={}; source_est={}
    treatments=('CT','L','F','L+F')
    for endpoint in ('drift','benthos'):
        for measure in ('abund','rich'):
            entries=outcomes[endpoint][measure]
            if len(entries)!=2 or entries[0]['treatment']!='L' or entries[1]['treatment']!='L+F':
                return {'status':'SOURCE_STOP_UNEXPECTED_CONTRAST_PAIR'}
            for entry in entries:
                for arm,field in ((entry['treatment'],'block_treatment_values'),
                                  (entry['control'],'block_control_values')):
                    key=f'{endpoint}|{measure}|{arm}'
                    vals=entry[field]
                    if len(vals)!=6 or any(not isinstance(x,(int,float)) or x<=0 or not math.isfinite(x) for x in vals):
                        return {'status':'SOURCE_STOP_INCOMPLETE_HT_COVARIANCE'}
                    if key in blockcolumns and blockcolumns[key]!=vals:
                        return {'status':'SOURCE_STOP_INCONSISTENT_SHARED_ARM'}
                    blockcolumns[key]=vals
                    means[key]=statistics.mean(vals)
            for idx,entry in enumerate(entries):
                ekey=f'{endpoint}|{measure}|{"L_vs_CT" if idx==0 else "LF_vs_F"}'
                names.append(ekey)
                source_est[ekey]=entry['lnRR']
                gradients.append({f'{endpoint}|{measure}|{entry["treatment"]}':1.0/entry['mean_treatment'],
                                  f'{endpoint}|{measure}|{entry["control"]}':-1.0/entry['mean_control']})
    # One influence-score vector per effect, with each original header equally
    # weighted. This generates a covariance of rank <= n_HT-1 even for 8 outputs.
    infl=[]
    for g in gradients:
        vals=[]
        for h in range(6):
            vals.append(sum(weight*(blockcolumns[key][h]-means[key])
                            for key,weight in g.items()))
        infl.append(vals)
    matrix=[[sum(infl[i][k]*infl[j][k] for k in range(6))/(6*5)
             for j in range(len(names))] for i in range(len(names))]
    tcrit=2.570581835636305
    def diff(a,b,kind):
        ai=names.index(a);bi=names.index(b)
        est=source_est[a]-source_est[b]
        var=matrix[ai][ai]+matrix[bi][bi]-2*matrix[ai][bi]
        if var<-1e-10:return {'status':'SOURCE_STOP_NEGATIVE_JOINT_VARIANCE'}
        se=math.sqrt(max(0.0,var))
        return {'kind':kind,'numerator_lnRR':a,'denominator_lnRR':b,
                'difference_lnRR':est,'ratio_of_response_ratios':math.exp(est),
                'SE_shared_HT':se,'t5_exploratory_95ratio':
                [math.exp(est-tcrit*se),math.exp(est+tcrit*se)],
                'independent_header_blocks':6,
                'status':'POSTHOC_DESCRIPTIVE_NO_CAUSAL_MEDIATION'}
    crosses=[]
    for metric in ('abund','rich'):
        for exposure in ('L_vs_CT','LF_vs_F'):
            crosses.append(diff(f'drift|{metric}|{exposure}',f'benthos|{metric}|{exposure}',
                                'DRIFT_VS_BENTHOS_COMPARISON_NOT_EMERGENCE'))
    flow_mod=[]
    for endpoint in ('drift','benthos'):
        for metric in ('abund','rich'):
            flow_mod.append(diff(f'{endpoint}|{metric}|LF_vs_F',f'{endpoint}|{metric}|L_vs_CT',
                                 'FLOW_MODIFIES_LIGHT_ASSOCIATION'))
    return {'status':'SOURCE_DERIVED_EIGHT_JOINT_HT_CONTRAST_COVARIANCE',
            'physical_programme':'EERES2021_Disconnected18flumes_6headers',
            'source_sha256':SOURCE_SHA,'independent_header_blocks':6,
            'source_outcomes':names,'covariance_logRR':matrix,
            'stage_contrasts':crosses,'light_by_flow_modifiers':flow_mod,
            'interpretation_gate':'Drift=in-stream downstream transport, benthos=riverbed community; neither is adult emergence or predation. Correlated outcomes across six source headers, 8 selected comparisons and 8 post hoc differences. No wet-vs-dry gradient or plant outcome.'}

def find_archive(archive):
    if len(archive)>MAX:raise ValueError('ZIP_TOO_BIG')
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        entries=z.infolist()
        if len(entries)>30:raise ValueError('ZIP_ENTRIES_TOO_MANY')
        candidates=[x for x in entries if x.filename.split('/')[-1]==NAME]
        if len(candidates)!=1:raise ValueError('SOURCE_WORKBOOK_NOT_UNIQUE')
        info=candidates[0]
        if info.file_size>1_000_000:raise ValueError('XLSX_TOO_BIG')
        raw=z.read(info)
        if raw[:2]!=b'PK':raise ValueError('NOT_REAL_XLSX')
        return {'filename':info.filename,'sha256':hashlib.sha256(raw).hexdigest(),
                'size_bytes':len(raw),'metadata':inspect_xlsx(raw),
                'source_effects':source_effects(raw) if hashlib.sha256(raw).hexdigest()==SOURCE_SHA else {'status':'SOURCE_STOP_HASH_CHANGED'}}

def run():
    result={'review':'EERES2021_Disconnected','doi':DOI,'journal_publications':[
     '10.1016/j.scitotenv.2023.168836','10.1002/ecog.07758'],
     'expected':NAME,'raw_record_republication':False,'effect_adopted':False,
     'status':'SOURCE_STOP','metadata_retrieved':False,'probes':[]}
    for endpoint in [B+'/datasets/doi:'+DOI,B+'/datasets/doi:'+quote(DOI,safe='')]:
        try:
            raw=grab(endpoint);obj=json.loads(raw)
            result['dataset_metadata_sha256']=hashlib.sha256(raw).hexdigest()
            result['dataset_title']=obj.get('title')
            result['dataset_id']=obj.get('id')
            result['metadata_retrieved']=True
            result['probes'].append({'url':endpoint,'status':'JSON_RETRIEVED'})
            break
        except Exception as e:result['probes'].append({'url':endpoint,'status':'SOURCE_STOP','error':type(e).__name__+':'+str(e)[:100]})
    else:return result
    targets=sorted(hrefs(obj));result['metadata_links']=targets[:20]
    version_links=sorted(u for u in targets if re.search(r'/versions/\d+$',urlparse(u).path))
    # Verify version is advertised by official dataset metadata, not guessed.
    for version in version_links:
        try:
            content=grab(version+'/download')
            result['zip_bytes']=len(content);result['zip_sha256']=hashlib.sha256(content).hexdigest()
            result['original_workbook']=find_archive(content)
            result['status']='SOURCE_WORKBOOK_VERIFIED'
            return result
        except Exception as e:result['probes'].append({'url':version+'/download','status':'SOURCE_STOP','error':type(e).__name__+':'+str(e)[:130]})
    result['status']='SOURCE_STOP_ARCHIVE';return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='results/eeres2021_disconnected_source_audit.json');args=ap.parse_args()
    d=run();p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('EERES_STATUS',d['status'])
    if 'original_workbook' in d:
        f=d['original_workbook'];print('SOURCE_WORKBOOK',f['sha256'],f['size_bytes'])
        for nm,s in f['metadata'].get('sheets',{}).items():
            print('SHEET',nm,s.get('n_nonempty_data_rows'),s.get('n_columns'))
    print('AUDIT',p)
