#!/usr/bin/env python3
"""Original Mendeley v1 ALAN flowering response at true 2023 site grain.

Grognuz et al., Biological Conservation 2026 doi:10.1016/j.biocon.2025.111650
Source: doi:10.17632/xjmn5jhxyn.1 (NOT presumed identical to reserved v2).
Report descriptive site-equal timing shifts and site-level resampling only,
not a reproduction of the final mixed model nor observed pollinator mismatch.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

SHA_DIEL="887d6f4b2a22216a9f08988c0eb30f03ed6a3d3476c1e5f6e7192a06579505a0"
SHA_SEASON="30b857acb6538a08298c039118cd54fd65e669af41a2c5bfcf4de75547c29f15"
EVENTS=[
  ("Silene noctiflora","Starts.closing","Time.of.sunrise"),
  ("Malva moschata","Starts.opening","Time.of.sunrise"),
  ("Silene noctiflora","Starts.opening","Time.of.sunset"),
  ("Malva moschata","Starts.closing","Time.of.sunset")
]

def check_source(path,expect):
    raw=Path(path).read_bytes()
    actual=hashlib.sha256(raw).hexdigest()
    if actual!=expect:raise ValueError("Original Mendeley v1 source SHA drift: "+actual)
    return raw

def clock_offset(time,reference):
    """Signed event-time shift relative to daily sunrise/sunset in minutes."""
    return ((time-reference+720)%1440)-720

def site_diff_bootstrap(lit,dark,seed,iterations=10000):
    """Site-means and independent SITE resampling, not original paired design."""
    x=np.asarray(lit,float);y=np.asarray(dark,float)
    if len(x)<2 or len(y)<2:return None
    rng=np.random.default_rng(seed)
    dx=rng.choice(x,size=(iterations,len(x)),replace=True).mean(axis=1)
    dy=rng.choice(y,size=(iterations,len(y)),replace=True).mean(axis=1)
    return np.quantile(dx-dy,[.025,.975]).tolist()

def diel_analysis(raw,iterations=10000):
    from io import BytesIO
    d=pd.read_csv(BytesIO(raw),sep=";")
    assert len(d)==555
    assert set(d.Species)=={"Silene noctiflora","Malva moschata"}
    site_names=d.groupby("Study.site").LIGHT_TRMT.nunique()
    assert len(site_names)==12 and site_names.eq(1).all()
    assert (d.drop_duplicates("Study.site").LIGHT_TRMT.value_counts().to_dict()=={"DARK":6,"LIGHT":6})
    out=[]
    for idx,(species,stage,sun) in enumerate(EVENTS):
        q=d[d.Species.eq(species)&d.Exp_Night.between(1,3)].copy()
        evt=pd.to_timedelta(q[stage],errors="coerce").dt.total_seconds()/60
        ref=pd.to_timedelta(q[sun],errors="coerce").dt.total_seconds()/60
        q["offset_min"]=clock_offset(evt,ref)
        q=q[q.offset_min.notna()]
        site=q.groupby(["Study.site","LIGHT_TRMT"],as_index=False).offset_min.median()
        lit=site[site.LIGHT_TRMT.eq("LIGHT")].offset_min.to_numpy(float)
        dark=site[site.LIGHT_TRMT.eq("DARK")].offset_min.to_numpy(float)
        assert (len(lit),len(dark))==(6,6)
        wk_site=q.groupby(["Week","Study.site","LIGHT_TRMT"],as_index=False).offset_min.median()
        weekly=wk_site.groupby(["Week","LIGHT_TRMT"]).offset_min.mean().unstack().dropna()
        out.append({"species":species,"event":stage,"reference":sun,
          "event_records":int(len(q)),"unique_plant_labels":int(q.Label.nunique()),
          "n_lit_sites":len(lit),"n_dark_sites":len(dark),
          "mean_of_site_medians_lit_minus_dark_minutes":float(lit.mean()-dark.mean()),
          "site_bootstrap95_descriptive_minutes":site_diff_bootstrap(lit,dark,20261115+idx,iterations),
          "n_matched_calendar_weeks":int(len(weekly)),
          "week_lit_later":int((weekly.LIGHT>weekly.DARK).sum()),
          "week_lit_earlier":int((weekly.LIGHT<weekly.DARK).sum()),
          "not_all_wild_Malva_data_in_original_paper":True,
          "week_and_plant_repeats_not_independent_site_replications":True})
    assert [(x["event_records"],x["unique_plant_labels"]) for x in out][:2]==[(137,68),(82,35)]
    assert abs(out[0]["mean_of_site_medians_lit_minus_dark_minutes"]-77.0833333333333)<1e-8
    assert abs(out[1]["mean_of_site_medians_lit_minus_dark_minutes"]+24.5)<1e-8
    return {"original_2023_potted_flower_day_rows":555,"n_true_light_assignment_sites":12,
        "short_term_source_only_not_all_91_wild_long_term_plants":True,
        "event_specific_results":out}

def seasonal_analysis(raw,iterations=10000):
    from io import BytesIO
    d=pd.read_csv(BytesIO(raw),low_memory=False)
    assert d.shape==(379080,18)
    assert (d.PLOT_ID.nunique(),d.PLANT_SP.nunique(),d.DOY.nunique(),d.YEAR_STUDY_SITE.nunique())==(108,27,130,27)
    fixed=d.groupby(["PLOT_ID","PLANT_SP"]).KEEP.nunique()
    assert len(fixed)==2916 and fixed.eq(1).all()
    flags=d.drop_duplicates(["PLOT_ID","PLANT_SP"]).KEEP.value_counts().to_dict()
    assert flags=={0:2066,1:850}
    assert d[d.FLORAL_UNIT.gt(0)].KEEP.eq(1).all()
    balance=d.drop_duplicates("YEAR_STUDY_SITE").groupby(["YEAR","LIGHT_TREATMENT"]).size().to_dict()
    assert balance[(2023,"LIT")]==8 and balance[(2023,"DARK")]==8
    assert balance[(2022,"LIT")]==6 and balance[(2022,"DARK")]==5
    observed=d[d.YEAR.eq(2023)&d.KEEP.eq(1)&d.FLORAL_UNIT.notna()]
    pos=observed[observed.FLORAL_UNIT.gt(0)]
    keys=["YEAR_STUDY_SITE","LIGHT_TREATMENT","PLANT_SP"]
    sample=pos.groupby(keys,as_index=False).DOY.agg(first="min",last="max")
    # A site is the experimental unit; aggregate the site's 1-4 plot
    # flower counts on each date before choosing that species' peak.
    daily=pos.groupby(keys+["DOY"],as_index=False).FLORAL_UNIT.sum()
    peak=(daily.sort_values(keys+["FLORAL_UNIT","DOY"],ascending=[True,True,True,False,True])
          .drop_duplicates(keys).loc[:,keys+["DOY"]].rename(columns={"DOY":"peak"}))
    sample=sample.merge(peak,on=keys,how="left",validate="one_to_one")
    first_obs=observed.groupby(keys,as_index=False).DOY.min().rename(columns={"DOY":"observed_first"})
    sample=sample.merge(first_obs,on=keys,validate="one_to_one")
    sample["first_record_left_censored"]=sample["first"].eq(sample["observed_first"])
    assert len(sample)==206
    support=[]
    for i,(sp,g) in enumerate(sorted(sample.groupby("PLANT_SP"),key=lambda t:t[0])):
        a=g[g.LIGHT_TREATMENT.eq("LIT")];b=g[g.LIGHT_TREATMENT.eq("DARK")]
        if not len(a) or not len(b):continue
        n_lit=len(a);n_dark=len(b)
        q={"species":sp,"n_lit_sites_flowering":int(n_lit),"n_dark_sites_flowering":int(n_dark),
           "n_lit_onset_left_censored":int(a.first_record_left_censored.sum()),
           "n_dark_onset_left_censored":int(b.first_record_left_censored.sum()),
           "site_support_ge4":bool(min(n_lit,n_dark)>=4)}
        for j,(src,name) in enumerate([("first","first"),("peak","peak"),("last","last")]):
            x=a[src].to_numpy(float);y=b[src].to_numpy(float)
            q[name+"_lit_minus_dark_days"]=float(x.mean()-y.mean())
            q[name+"_site_resample95"]=site_diff_bootstrap(x,y,20261115+100*i+j,iterations) if q["site_support_ge4"] else None
        support.append(q)
    assert sum(x["site_support_ge4"] for x in support)==12
    sensitivity={}
    for threshold in (3,4,5,6):
        q=[x for x in support if x["n_lit_sites_flowering"]>=threshold and x["n_dark_sites_flowering"]>=threshold]
        sensitivity[str(threshold)]={"species":len(q),
          "onset_earlier_n":sum(x["first_lit_minus_dark_days"]<0 for x in q),
          "peak_earlier_n":sum(x["peak_lit_minus_dark_days"]<0 for x in q)}
    assert sensitivity["4"]=={"species":12,"onset_earlier_n":7,"peak_earlier_n":10}
    assert sensitivity["5"]=={"species":6,"onset_earlier_n":5,"peak_earlier_n":5}
    return {
      "original_2022_2023_source_plot_species_day_expanded_rows":int(len(d)),
      "original_unique_plot_species_cohorts":int(len(fixed)),
      "source_keep1_plot_species_cohorts":850,
      "source_keep0_plot_species_cohorts":2066,
      "positive_flower_count_row_n":int(d.FLORAL_UNIT.gt(0).sum()),
      "year_site_counts":{f"{yr}_{light}":int(count) for (yr,light),count in balance.items()},
      "analyzed_year":2023,"n_source_2023_site_species_positive":int(len(sample)),
      "onset_positive_at_first_observed_n":int(sample.first_record_left_censored.sum()),
      "onset_positive_at_first_observed_by_light":{
        light:int(g.first_record_left_censored.sum()) for light,g in sample.groupby("LIGHT_TREATMENT")},
      "site_support_threshold_sensitivities":sensitivity,
      "species_site_effects":support,
      "flowering_cohort_conditioned_on_observed_positive_flower_counts":True,
      "missing_records_not_zero_flower_reclassified":True,
      "reconstructed_first_peak_last_are_not_authors_full_mixed_model":True,
      "species_level_CIs_not_multiplicity_adjusted":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--diel",required=True)
    ap.add_argument("--seasonal",required=True)
    ap.add_argument("--out",required=True)
    ap.add_argument("--iterations",type=int,default=10000)
    args=ap.parse_args()
    dst=Path(args.out)
    if dst.exists():raise FileExistsError("Refuse to overwrite source result")
    original_diel=check_source(args.diel,SHA_DIEL)
    original_season=check_source(args.seasonal,SHA_SEASON)
    result={"schema":"ALAN_Grognuz2026_source_v1_site_equal_diel_seasonal_v1",
      "paper_DOI":"10.1016/j.biocon.2025.111650",
      "public_source_DOI":"10.17632/xjmn5jhxyn.1",
      "paper_reserved_v2_NOT_verified_identical":True,
      "source_SHA256":{"diel":SHA_DIEL,"seasonal":SHA_SEASON},
      "site_pair_identity_not_assigned_from_animal_name_codes":True,
      "no_measured_pollinator_daily_or_seasonal_overlap_in_source":True,
      "no_seed_fitness_outcome_in_source":True,
      "diel":diel_analysis(original_diel,args.iterations),
      "seasonal":seasonal_analysis(original_season,args.iterations)}
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GROGNUZ2026_TRUE_SITE_BASED_BIOLOGY",json.dumps({
      "diel":[{"sp":v["species"],"event":v["event"],"n":v["event_records"],"sites":v["n_lit_sites"]+v["n_dark_sites"],"delta_min":v["mean_of_site_medians_lit_minus_dark_minutes"],"CI95":v["site_bootstrap95_descriptive_minutes"]} for v in result["diel"]["event_specific_results"]],
      "seasonal_thresholds":result["seasonal"]["site_support_threshold_sensitivities"],
      "source_2023_site_species":result["seasonal"]["n_source_2023_site_species_positive"]
      },ensure_ascii=False))
if __name__=="__main__":main()
