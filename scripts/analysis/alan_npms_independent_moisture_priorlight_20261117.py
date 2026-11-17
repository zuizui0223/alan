#!/usr/bin/env python3
"""Independent UK NPMS repeated *positive recorded cover* × prior PCNL × species MOIST.

Source: original signed PR#31 artifact 13763979354 + official Tichy EIV Europe
CSV, verified with SHA. Analysis uses plants regardless of ALAN-283 candidate status.
No zero is inferred for unreported plants; no causal claim or site-year exposure to lux.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr
from scipy.stats import t

EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
SOURCE_ORIGINAL_PRIOR_PCNLSHA="5ecae6a9b62de636c8e68fd694c69698eec154e5a8a0ba641beb5a7e6061e427"
SCHEMA="alan_npms_moisture_alan_independent_plant_cover_v1"
MIN_ORIG_RANGE=.02
MIN_SUPPORT={"rows":2500,"plots":150,"species":40,"wet_species":10,
             "drier_species":20,"mixed_plotyears":150}

def assert_source(ok,msg):
    if not bool(ok): raise ValueError("ALAN_NPMS_MOISTURE_SOURCE: "+str(msg))

def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for q in iter(lambda:f.read(4<<20),b""):h.update(q)
    return h.hexdigest()

def load(original, receipt, eiv):
    audit=json.loads(Path(receipt).read_text())
    assert_source(audit["source_original_positive_cover_table_SHA256"]==sha256(original),
                  "original author-pinned positive plant cover bytes changed")
    assert_source(audit.get("original_16960_plotyear_light_SHA256")==SOURCE_ORIGINAL_PRIOR_PCNLSHA,
                  "2014-2024 original PCNL plotyear SHA does not match")
    assert_source(audit.get("source_CRC_pinned_DOMIN_midPoint_unit")==
                  "fraction of cover [0,1], not percent [0,100]","DOMIN lookup basis changed")
    assert_source(audit.get("source_2015to2025_same_surveyID_season_eligible_unique_plots")==1203,
                  "original NPMS eligible plot frame changed")
    assert_source(sha256(eiv)==EIV_SHA,"official Tichy EIV source SHA changed")
    z=pd.read_csv(original,low_memory=False,dtype={"location_id":str,"survey_method_id":str,"species":str,"season":str})
    req={"location_id","survey_method_id","season","species","year",
         "median_cover","median_source_day_of_year","source_domin_encodings",
         "source_record_status","original_annual_priorPCNL","source_prior_unique_values"}
    assert_source(req.issubset(z.columns),"original plot×year×plant columns absent")
    panel_cols=["location_id","survey_method_id","season","species","year"]
    assert_source(not z.duplicated(panel_cols).any(),"pseudo-replicated source panel")
    assert_source(z.year.between(2015,2025).all() and
                  z.median_cover.between(0,1).all() and
                  z.median_source_day_of_year.between(1,366).all() and
                  z.original_annual_priorPCNL.ge(0).all() and
                  z.source_prior_unique_values.eq(1).all(),"invalid source outcome/light/year")
    z["panel_key"]=z[panel_cols[:4]].agg("|".join,axis=1)
    z["prior_light"]=z.original_annual_priorPCNL.astype(float)
    grp=z.groupby("panel_key",sort=False).agg(
        yr=("year","nunique"),lvals=("prior_light","nunique"),
        lo=("prior_light","min"),hi=("prior_light","max"),
        domin=("source_domin_encodings","nunique"))
    allowed=grp.loc[(grp.yr>=2)&(grp.lvals>=2)&((grp.hi-grp.lo)>=MIN_ORIG_RANGE)&grp.domin.eq(1)].index
    z=z[z.panel_key.isin(allowed)].copy()
    assert_source((len(z),z.location_id.nunique(),z.species.nunique())==(61555,1138,686),
                  "previously replicated positive cover source cohort changed")
    # Same published Europe EIV table as current main. No fuzzy name replacement.
    ei=pd.read_csv(eiv,encoding="latin1",index_col=0,low_memory=False)
    assert_source("MOIST" in ei.columns,"moisture column not in original indicator table")
    ei["MOIST"]=pd.to_numeric(ei["MOIST"],errors="coerce")
    # If the original source has conflicting same-name rows, do not choose one.
    dupl=ei.groupby(level=0).MOIST.agg(["first","nunique"])
    m=dupl.loc[dupl["nunique"]<=1,"first"].dropna()
    z["moist"]=z.species.map(m)
    # The lookup is exact published binomial; synonyms remain missing, not guessed.
    source_info={"original_source_rows":len(z),
      "original_source_plots":int(z.location_id.nunique()),
      "original_source_species":int(z.species.nunique()),
      "moisture_source_exact_taxa":int(z.loc[z.moist.notna(),"species"].nunique()),
      "moisture_source_exact_rows":int(z.moist.notna().sum()),
      "eiv_original_sha256":sha256(eiv),
      "source_original_cover_sha256":sha256(original)}
    z=z[z.moist.between(1,12)].copy()
    z["wet7"]=z.moist.gt(7).astype(float)
    z["moist_centered"]=z.moist.astype(float)-7.
    z["light_moist"]=z.prior_light*z.moist_centered
    z["light_wet7"]=z.prior_light*z.wet7
    z["doy30"]=(z.median_source_day_of_year.astype(float)-180.)/30.
    z["cover_pp"]=100*z.median_cover.astype(float)
    z["species_year"]=z.species+"|"+z.year.astype(str)
    z["plotyear_method_season"]=(
        z.location_id+"|"+z.survey_method_id+"|"+z.season+"|"+z.year.astype(str))
    assert_source(np.isfinite(z[["moist","prior_light","doy30","cover_pp"]]).all().all(),
                  "nonfinite known original measurements")
    return z,source_info

def eligibility(z):
    # Common support of distinct moisture guilds in the exact same
    # physical plot×survey×season×year that shares source satellite light.
    wide=z.groupby("plotyear_method_season").agg(
        wet=("wet7","max"),dry=("wet7","min"),
        rows=("species","nunique"),plots=("location_id","first"),years=("year","first"))
    mixed=wide[(wide.wet==1)&(wide.dry==0)&(wide.rows>=2)]
    out={"source_scored_rows":len(z),"source_scored_plots":int(z.location_id.nunique()),
         "source_scored_species":int(z.species.nunique()),
         "wet_species":int(z.loc[z.wet7==1,"species"].nunique()),
         "drier_species":int(z.loc[z.wet7==0,"species"].nunique()),
         "mixed_source_plotyear_methodseason_units":len(mixed),
         "mixed_source_plots":int(mixed.plots.nunique()),
         "mixed_source_calendar_years":int(mixed.years.nunique()),
         "original_source_same_2015_2025_method_season_light_span_min":MIN_ORIG_RANGE}
    out["pass_primary"]=(out["source_scored_rows"]>=MIN_SUPPORT["rows"] and
        out["source_scored_plots"]>=MIN_SUPPORT["plots"] and
        out["source_scored_species"]>=MIN_SUPPORT["species"] and
        out["wet_species"]>=MIN_SUPPORT["wet_species"] and
        out["drier_species"]>=MIN_SUPPORT["drier_species"])
    out["pass_strict_same_plotyear"]=(out["pass_primary"] and
        out["mixed_source_plotyear_methodseason_units"]>=MIN_SUPPORT["mixed_plotyears"])
    return mixed,out

def within_residualize(z,y,terms,fe_groups,w,eps=2e-7):
    # Residualizes all columns in the weighted metric against original
    # plot-plant and species-year, and optionally source plot-survey-year FEs.
    n=len(z)
    arrays=[pd.factorize(z[f].to_numpy(),sort=True)[0].astype(np.int32) for f in fe_groups]
    ns=[int(x.max()+1) for x in arrays]
    rr=np.tile(np.arange(n,dtype=np.int32),len(arrays))
    offsets=np.cumsum([0]+ns[:-1])
    cc=np.concatenate([a+o for a,o in zip(arrays,offsets)])
    F=csr_matrix((np.ones(len(rr),dtype=float),(rr,cc)),shape=(n,sum(ns)))
    s=np.sqrt(w)
    FW=F.multiply(s[:,None]).tocsr()
    Z=z[[y,*terms]].to_numpy(float)
    out=np.empty_like(Z)
    diags=[]
    for i in range(Z.shape[1]):
        res=lsmr(FW,Z[:,i]*s,atol=eps,btol=eps,maxiter=3300)
        assert_source(res[1] not in (0,3,6,7) and res[2]<3300,
                      f"source FE nonconvergence column={i} istop={res[1]} it={res[2]}")
        out[:,i]=Z[:,i]-F@res[0]
        diags.append({"col":i,"iterations":int(res[2]),"termination":int(res[1]),
                      "weighted_residual_l2":float(np.linalg.norm(out[:,i]*s))})
    return out,{"levels":dict(zip(fe_groups,ns)),"solver":diags}

def fit(z,terms,fe_names,weight_mode="plot_equal"):
    z=z.copy()
    n=len(z);nplot=z.location_id.nunique()
    assert_source(n>500 and nplot>=35, "too few original observed plants/plots for model")
    if weight_mode=="plot_equal":
        w=1./z.groupby("location_id").location_id.transform("size").to_numpy(float)
    elif weight_mode=="record_equal":w=np.ones(n)
    else:raise ValueError("Unknown weighting")
    v,fe=within_residualize(z,"cover_pp",terms,fe_names,w)
    Y,X=v[:,0],v[:,1:]
    S=(X*w[:,None]).T@X
    cond=float(np.linalg.cond(S))
    if not np.isfinite(cond) or cond>1e10 or np.linalg.matrix_rank(S)!=len(terms):
        return {"status":"no_identifiable_term_after_true_source_FE","rows":n,
                "plots":int(nplot),"condition_number":cond,"fe":fe}
    B=np.linalg.inv(S); beta=B@((X*w[:,None]).T@Y)
    u=Y-X@beta
    clusters=pd.factorize(z.location_id)[0]
    score=np.zeros((nplot,len(terms)),float)
    np.add.at(score,clusters,X*(w*u)[:,None])
    V=B@(score.T@score)@B
    V*=nplot/(nplot-1)*n/(n-len(terms))
    se=np.sqrt(np.maximum(np.diag(V),0))
    if not np.isfinite(se).all() or (se<=0).any():
        return {"status":"plot_cluster_SE_unidentified","rows":n,"plots":int(nplot),"fe":fe}
    crt=t.ppf(.975,nplot-1)
    return {"status":"observational_positive_cover_only", "rows":n,"plots":int(nplot),
        "taxa":int(z.species.nunique()),"original_years":int(z.year.nunique()),
        "FE":fe_names,"weight":weight_mode,"condition_number":cond,"fe":fe,
        "coefficients":{key:{"beta":float(beta[i]),"SE_real_plot_cluster":float(se[i]),
             "nominal_95CI":[float(beta[i]-crt*se[i]),float(beta[i]+crt*se[i])],
             "nominal_two_sided_p":float(2*t.sf(abs(beta[i]/se[i]),nplot-1))}
             for i,key in enumerate(terms)}}

def execute(original,receipt,eiv):
    z,info=load(original,receipt,eiv)
    mixed,support=eligibility(z)
    result={"schema":SCHEMA,"source":info,"eligibility":support,
            "interpretation":"independent NPMS 2015-25 positive-only recorded plant cover; exact Europe EIV match",
            "endpoint":"DOMIN measured positive cover percentage points; non-records not coded zero",
            "hypothesis":"Prior annual satellite PCNL x independent species MOIST moderation",
            "coefficients_not_causal":True,"not_plant_loss_or_extinction":True,
            "frozen_thresholds":MIN_SUPPORT}
    if not support["pass_primary"]:
        result["stop_reason"]="insufficient_source_species_moisture_and_plot_support"
        return result
    # Primary direct test: all scored positively recorded plants, with same
    # physical plot-plant-method-season and true plant x year fixed effects.
    terms=["prior_light","light_moist","doy30"]
    F=["panel_key","species_year"]
    result["primary_continuous_moisture"]=fit(z,terms,F)
    # Explicit categorical no-cutoff-search sensitivity, only MOIST >7
    result["wet7_sensitivity"]=fit(z,["prior_light","light_wet7","doy30"],F)
    if not support["pass_strict_same_plotyear"]:
        result["strict_stop_reason"]="not enough truly co-observed wet and drier taxa"
        return result
    # Strict comparison controls *each identical plot×method×season×year*.
    # Common yearly light level is absorbed, but interaction may vary across
    # co-occurring wet and dry species with different observational trajectories.
    keep=z[z.plotyear_method_season.isin(mixed.index)].copy()
    result["strict_source_rows"]=len(keep)
    result["strict_moisture"]=fit(
         keep,["light_moist","doy30"],
         ["panel_key","species_year","plotyear_method_season"])
    result["strict_wet7"]=fit(
         keep,["light_wet7","doy30"],
         ["panel_key","species_year","plotyear_method_season"])
    return result

def cli():
    a=argparse.ArgumentParser()
    a.add_argument("--source",required=True)
    a.add_argument("--receipt",required=True)
    a.add_argument("--eiv",required=True)
    a.add_argument("--out",required=True)
    args=a.parse_args()
    result=execute(args.source,args.receipt,args.eiv)
    p=Path(args.out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="source"},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    cli()
