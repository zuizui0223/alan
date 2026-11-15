"""Actual same-site/date inter-flower timing window, not pollinator mismatch.

Two true taxa at 12 independently illuminated/dark sites (6+6).
Coavailability is BOUNDED by flower movement onset versus fully-open/closed.
Site-date may contain different individual plants, so this is a modelled
pair-of-species window, NOT confirmed simultaneous effective pollen receipt.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_SHA="887d6f4b2a22216a9f08988c0eb30f03ed6a3d3476c1e5f6e7192a06579505a0"
STAGES=[
 ("Silene noctiflora","Starts.closing","silene_begin_close"),
 ("Silene noctiflora","Fully.closed","silene_full_closed"),
 ("Malva moschata","Starts.opening","malva_begin_open"),
 ("Malva moschata","Fully.open","malva_full_open")]
BOUNDS={
 "permissive_partially_open":("silene_full_closed","malva_begin_open"),
 "both_transition_begun":("silene_begin_close","malva_begin_open"),
 "malva_fully_open_silene_not_closed":("silene_full_closed","malva_full_open"),
 "conservative_both_fully_open":("silene_begin_close","malva_full_open")}
WINDOWS=[90,120,150,180,240]

def minutes(s):
    return pd.to_timedelta(s,errors="coerce").dt.total_seconds()/60

def circadian_offset(t,reference):
    """Place clock time on nearest day relative to sunrise in [-12h,12h)."""
    return ((t-reference+720)%1440)-720

def dawn_window(close,op,half_width=150):
    """Modelled co-available minutes in [sunrise-half_width,sunrise+half_width]."""
    return np.maximum(0,np.minimum(close,half_width)-np.maximum(op,-half_width))

def site_effect(z,close,op,half_width,seed=20261115,B=15000):
    x=z.copy()
    x["minutes"]=dawn_window(x[close].to_numpy(float),x[op].to_numpy(float),half_width)
    site=x.groupby(["Study.site","LIGHT_TRMT"],as_index=False).minutes.mean()
    a=site.loc[site.LIGHT_TRMT.eq("LIGHT"),"minutes"].to_numpy(float)
    b=site.loc[site.LIGHT_TRMT.eq("DARK"),"minutes"].to_numpy(float)
    assert len(a)==len(b)==6
    rng=np.random.default_rng(seed)
    sample=(rng.choice(a,(B,6),replace=True).mean(axis=1)
            -rng.choice(b,(B,6),replace=True).mean(axis=1))
    return {
       "observed_site_days":int(len(x)),
       "lit_site_days":int(x.LIGHT_TRMT.eq("LIGHT").sum()),
       "dark_site_days":int(x.LIGHT_TRMT.eq("DARK").sum()),
       "independent_lit_sites":6,"independent_dark_sites":6,
       "lit_site_equal_mean_minutes":float(a.mean()),
       "dark_site_equal_mean_minutes":float(b.mean()),
       "lit_minus_dark_minutes":float(a.mean()-b.mean()),
       "site_resample95_descriptive_minutes":np.quantile(sample,[.025,.975]).tolist(),
       "n_lit_sites_with_any_overlap":int((a>0).sum()),
       "n_dark_sites_with_any_overlap":int((b>0).sum()),
       "stage_proxy_not_directly_observed_both_open":True,
       "site_labels_not_recovered_as_randomized_pairs":True
    }

def analyze(data):
    assert len(data)==555
    assert set(data.Species)=={"Silene noctiflora","Malva moschata"}
    q=data[data.Exp_Night.between(1,3)].copy()
    assert q["Study.site"].nunique()==12
    outs={};counts={}
    for species,source_field,key in STAGES:
        t=q[q.Species.eq(species)].copy()
        offset=circadian_offset(minutes(t[source_field]),minutes(t["Time.of.sunrise"]))
        t[key]=offset
        t=t[t[key].notna()]
        counts[key]={
          "event_rows":len(t),
          "early_relative_to_sunrise_minus6h":int(t[key].lt(-360).sum()),
          "late_relative_to_sunrise_plus6h":int(t[key].gt(360).sum()),
          "late_or_early_by_treatment":{
           treatment:int((t.loc[t.LIGHT_TRMT.eq(treatment),key].abs()>360).sum())
           for treatment in ["LIGHT","DARK"]}}
        outs[key]=t.groupby(["Study.site","LIGHT_TRMT","Date"],as_index=False)[key].median()
    key=["Study.site","LIGHT_TRMT","Date"]
    z=outs["silene_begin_close"]
    for stage in [x[2] for x in STAGES[1:]]:
        z=z.merge(outs[stage],on=key,validate="one_to_one")
    assert len(z)==46 and z["Study.site"].nunique()==12
    assert z.LIGHT_TRMT.value_counts().to_dict()=={"DARK":24,"LIGHT":22}
    assert counts["silene_begin_close"]["event_rows"]==137
    assert counts["malva_begin_open"]["event_rows"]==82
    out={}
    for i,(name,(cl,op)) in enumerate(BOUNDS.items()):
       out[name]={str(w):site_effect(z,cl,op,w,seed=20261115+i*13+w)
                  for w in WINDOWS}
    # A gate for analytical selection: discarding flowers that start closing
    # >6h before sunrise preferentially excludes dark-site early closures.
    filtered=z[z.silene_begin_close.between(-360,360)&
               z.malva_begin_open.between(-360,360)]
    assert len(filtered)==42 and filtered.LIGHT_TRMT.value_counts().to_dict()=={"DARK":21,"LIGHT":21}
    strict_selected=site_effect(filtered,"silene_begin_close","malva_begin_open",150)
    result={
       "schema":"ALAN_Grognuz2026_true_site_date_crossflower_window_v1",
       "original_paper_DOI":"10.1016/j.biocon.2025.111650",
       "data_DOI":"10.17632/xjmn5jhxyn.1",
       "original_author_source_SHA256":SOURCE_SHA,
       "input_source_event_rows":len(data),
       "actual_silene_close_event_count":counts["silene_begin_close"]["event_rows"],
       "actual_malva_open_event_count":counts["malva_begin_open"]["event_rows"],
       "matched_true_site_date_pairs":int(len(z)),
       "source_record_stage_counts":counts,
       "within_site_date_source_plant_individuals_not_the_same":True,
       "stage_bound_definitions":{
        k:{"silene_stopping_condition":v[0],"malva_start_condition":v[1]}
        for k,v in BOUNDS.items()},
       "dawn_windows_minutes_each_side":WINDOWS,
       "site_equal_stage_window_sensitivity":out,
       "early_closure_selection_sensitivity":{
         "selection_rule":"exclude events over six hours from sunrise (not an a priori biological exclusion)",
         "n_matched_site_dates":int(len(filtered)),
         "started_transitions_only_5hour_window":strict_selected,
         "original_dark_early_closures_excluded_more_frequently":True},
       "H3_true_plant_pollinator_timing_mismatch_estimable":False,
       "pollinator_hourly_activities_measured":False,
       "same_flower_actual_petals_open_measured_in_this_join":False,
       "meaning":"potential cross-species floral availability timing, NOT effective pollen/seed service",
       "no_preferred_effect_sign_claimed":True}
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    raw=a.source.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:
        raise ValueError("Mendeley version 1 source byte identity changed")
    result=analyze(pd.read_csv(io.BytesIO(raw),sep=";"))
    if a.out.exists():raise FileExistsError("no overwrite")
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    focus={k:v["150"] for k,v in result["site_equal_stage_window_sensitivity"].items()}
    print("ALAN_SOURCE_DAWN_BOUND_RESULTS",json.dumps({
        "site_dates":result["matched_true_site_date_pairs"],
        "bounds":focus,
        "strict_window_selection_delta":result["early_closure_selection_sensitivity"]["started_transitions_only_5hour_window"]["lit_minus_dark_minutes"],
        "pollinator_overlap_estimated":False},ensure_ascii=False))
if __name__=="__main__":main()
