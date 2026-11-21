#!/usr/bin/env python3
"""Original NPMS observation-process and management-label robustness of ALAN×MOIST×habitat.

This does NOT infer that a survey with an unlisted plant means local extinction.
NPMS species checklists are habitat-dependent; management/grazing labels are
surveyor-reported observations, not interventions. Stable-label restrictions
diagnose how much of the wetness×PCNL signal relies on changes in observation
protocol and reported management, not an unbiased causal management adjustment.

All volunteer sample/site IDs stay in RAM. JSON exports aggregate counts only.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.acquisition import alan_npms_original_habitat_source_gate_20261117 as h
from scripts.analysis import alan_npms_true_wetland_moisture_priorlight_20261117 as wet
from scripts.analysis import alan_npms_wetland_light_spatial_falsification_20261117 as fl

SCHEMA="alan_npms_true_observer_management_and_checklist_gate_v1"
CAPTIONS=("NPMS Habitat","Management","Grazing","Grazers","Other Management","No species observed")
CHECKS=("NPMS Habitat","Management","Grazing","Grazers","Other Management")
def fail(ok,msg):
    if not bool(ok):raise ValueError("ALAN_NPMS_OBSERVER_PROCESS_STOP: "+msg)
def top(df):
    return {str(k):int(v) for k,v in df.value_counts(dropna=False).head(16).items()}
def source_unit_metadata(orig,common):
    sam=orig["sample"]
    att=orig["attrs"]
    assert sam.id.is_unique
    df=sam[["id","survey_id","location_id","date_start"]].copy()
    dt=pd.to_datetime(df.date_start,errors="coerce",utc=True)
    df["year"]=dt.dt.year
    df["season"]=np.select([dt.dt.month.isin([4,5,6]),dt.dt.month.isin([7,8,9])],
                          ["spring","summer"],default="other")
    df["survey_method_id"]=df.survey_id.str.strip()
    df=df[df.location_id.isin(set(common.location_id))].copy()
    keys=["location_id","year","season","survey_method_id"]
    # One NPMS habitat or management note is a property of a *survey visit*.
    # Multiple visits within plot-year remain separate until labels are
    # reconciled. A missing caption is UNKNOWN, not stable "No".
    att=att[att.caption.isin(CAPTIONS)].copy()
    att=att.merge(df[["id",*keys]],left_on="sample_id",right_on="id",
                  how="inner",validate="many_to_one")
    def collapse(z):
        vals=sorted(set(v.strip() for v in z if v.strip()))
        return "|".join(vals) if vals else None
    snap=att.groupby([*keys,"caption"],dropna=False).text_value.agg(collapse).unstack("caption")
    # When two visits for same survey unit have distinct categorical labels,
    # retain both as explicit 'inconsistent' rather than choosing one.
    snap=snap.reset_index()
    source_keys=common[keys].drop_duplicates()
    unit=source_keys.merge(snap,on=keys,how="left",validate="one_to_one")
    fail(len(unit)==len(source_keys),"nonunique true survey units")
    return att,unit

def stable_plot_labels(unit,col):
    d=unit[["location_id",col]].copy()
    # Eligibility requires every included source survey year to contain
    # a nonempty observed original source label, with exactly one unique
    # label within the repeated series. No filling from past observations.
    d[col]=d[col].replace("",np.nan)
    group=d.groupby("location_id",sort=False)[col]
    summary=group.agg(n_survey_units="size",n_with_label="count",n_labels="nunique")
    good=summary[(summary.n_survey_units>=2)&
                 (summary.n_with_label==summary.n_survey_units)&
                 (summary.n_labels==1)]
    return set(good.index),{
        "source_physical_plots":int(len(summary)),
        "plots_all_years_labeled":int((summary.n_with_label==summary.n_survey_units).sum()),
        "plots_stable_all_years_labeled":int(len(good)),
        "plots_with_two_or_more_observed_labels":int((summary.n_labels>=2).sum()),
        "plots_with_any_missing_label":int((summary.n_with_label<summary.n_survey_units).sum()),
    }

def choose_fit(common,site_ids,label):
    d=common[common.location_id.isin(site_ids)].copy()
    d=fl.mixed_reselect(d)
    x=fl.estimate(d,label)
    out={"status":x["status"],"source_support":x["support"],"retained_site_n":len(site_ids)}
    if "triple_habitat_moisture_light" in x:
        y=x["triple_habitat_moisture_light"]
        out.update(beta=y["beta"],plot_cluster_CI=y["plot_95CI"],
                   two_way_plot_species=y["two_way"])
    return out

def investigate(npms,cover,receipt,eiv):
    summary,common,eligible=wet.study(npms,cover,receipt,eiv,return_source_panels=True)
    fail(summary["source_gate"]["pass"],"source-locked habitat gate changed")
    fail(len(common)==11902 and common.location_id.nunique()==326,
         "source-locked 47 wetland + 279 controls changed")
    orig=h.read(npms)
    att,unit=source_unit_metadata(orig,common)
    output={
        "schema":SCHEMA,
        "dataset_NPMS_source":"10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8",
        "source_original_internal_CRC_verified":True,
        "original_cover_sha256":h.sha(cover),
        "source_plant_positive_cover_rows":int(len(common)),
        "wetland_source_eligible_plots":47,"other_source_eligible_plots":279,
        "observed_survey_units":int(len(unit)),
        "source_aggregate_captions":{},
        "source_attribute_coverage":{},
        "original_list_changes_are_observation_not_absence":True,
        "no_biological_zero_imputation":True,
        "no_assigned_causal_light_or_management_effect":True,
        "robustness":{}
    }
    for c in CAPTIONS:
        s=unit[c] if c in unit else pd.Series(dtype=str)
        output["source_attribute_coverage"][c]={
            "source_units_with_observed_label":int(s.notna().sum()),
            "source_units_without_label":int(s.isna().sum())
        }
        # Original values are text tokens, not invented numerical intensity.
        output["source_aggregate_captions"][c]=(top(s) if len(s) else {}) if c!="Other Management" else {"potentially_free_text_values_not_exported":int(s.notna().sum())}
    for c in CHECKS:
        retained,stats=stable_plot_labels(unit,c)
        summary_fit=choose_fit(common,retained,"stable_original_"+c.replace(" ","_"))
        output["robustness"][c]={"stability":stats,"model":summary_fit,
            "note":"conditioned on original nonmissing identical reported labels across surveyed years"}
    # Independently documented candidate pathways from the source *reported*
    # management descriptions. A note marked "water regime regulation" may
    # indicate hydrology/sluice management, not treatment ALAN. Conversely,
    # missing notes and other management labels are NOT evidence of no water
    # management. Keep aggregate matched site-year counts, not text or IDs.
    membership=common[["location_id","inland_wetland"]].drop_duplicates()
    fail(membership.location_id.is_unique,
         "same source physical plot has inconsistent wetland category")
    details=unit.merge(membership,on="location_id",validate="many_to_one")
    flag_terms={
        "water_regime_reported":"Water regime regulation",
        "ditch_clearance_reported":"Ditch-clearance",
        "cutting_mowing_reported":"Cutting / mowing",
        "livestock_grazing_reported":"Grazing - livestock",
        "wild_herbivore_grazing_reported":"Grazing - rabbits",
        "fenced_grazers_excluded_reported":"Fenced to exclude grazing"
    }
    output["direct_source_management_hydrology_candidates"]={}
    for wetflag,grp in details.groupby("inland_wetland"):
        label="inland_wetland" if wetflag else "other_habitat"
        valid=grp[grp["Management"].notna()].copy()
        body={"physical_plots":int(grp.location_id.nunique()),
              "all_source_plotyear_units":len(grp),
              "survey_units_with_observed_management":len(valid),
              "survey_units_without_management":len(grp)-len(valid)}
        for key,term in flag_terms.items():
            q=valid["Management"].str.contains(term,regex=False,na=False)
            body[key]={
              "source_survey_units_with_reported_code":int(q.sum()),
              "source_distinct_plots_with_reported_code":int(
                  valid.loc[q,"location_id"].nunique()),
              "source_units_with_other_management_code_not_proof_of_absence":int((~q).sum()),
              "source_eligible_missing_management_labels":len(grp)-len(valid)
            }
        output["direct_source_management_hydrology_candidates"][label]=body
    # Falsify a light-specific story against independently observed original
    # hydrological/vegetation management. A physical plot is excluded if
    # **any** of its source plotyears explicitly reported that management
    # code. No report != no management; these are incomplete-notes
    # sensitivity exclusions, NOT an adjustment for true water status.
    source_site_ids=set(common.location_id)
    valid_notes=details[details["Management"].notna()].copy()
    exclusions={
        "water_regulation":("Water regime regulation",),
        "ditch_clearance":("Ditch-clearance",),
        "either_water_regulation_or_ditch_clearance":
            ("Water regime regulation","Ditch-clearance"),
        "cutting_or_mowing":("Cutting / mowing",),
        "livestock_grazing":("Grazing - livestock",)
    }
    output["source_reported_management_exclusion_robustness"]={}
    for label,names in exclusions.items():
        source_flag=pd.Series(False,index=valid_notes.index)
        for term in names:
            source_flag |= valid_notes["Management"].str.contains(
                term,regex=False,na=False)
        flagged=set(valid_notes.loc[source_flag,"location_id"])
        checked=source_site_ids-flagged
        wet_sites=set(membership.loc[membership.inland_wetland,"location_id"])
        other_sites=source_site_ids-wet_sites
        fitted=choose_fit(common,checked,"exclude_source_reported_"+label)
        output["source_reported_management_exclusion_robustness"][label]={
            "source_reported_management_codes":list(names),
            "flagged_total_plots":len(flagged),
            "flagged_inland_wetland_plots":len(flagged&wet_sites),
            "flagged_other_plots":len(flagged&other_sites),
            "other_plots_with_missing_notes_are_kept_and_NOT_called_unmanaged":True,
            "status":fitted["status"],
            "source_support":fitted["source_support"],
            "beta":fitted.get("beta"),
            "nominal_95CI":fitted.get("plot_cluster_CI"),
            "two_way":fitted.get("two_way_plot_species")
        }
    # Exclude original plot x plant x method x season panels with changes in
    # NPMS record_status (C/V), which are not guaranteed comparable record
    # processes. A status change cannot be overwritten or called plant death.
    p=common.groupby("panel_key").source_record_status.nunique(dropna=False)
    stable_keys=set(p.index[p.eq(1)])
    status_data=common[common.panel_key.isin(stable_keys)].copy()
    status_data=fl.mixed_reselect(status_data)
    st_fit=fl.estimate(status_data,"stable_original_plant_record_status")
    output["robustness"]["stable_original_plant_record_status"]={
        "original_panel_count":int(p.size),
        "panels_stable_status":int(p.eq(1).sum()),
        "panels_with_status_switch":int(p.gt(1).sum()),
        "source_support":st_fit["support"],"status":st_fit["status"],
        "beta":st_fit.get("triple_habitat_moisture_light",{}).get("beta"),
        "nominal_CI":st_fit.get("triple_habitat_moisture_light",{}).get("plot_95CI"),
        "two_way":st_fit.get("triple_habitat_moisture_light",{}).get("two_way")
    }
    # Union of original nonmissing Management and Grazing stability is
    # deliberately strict and may fail support. Do not relax based on sign.
    a,_=stable_plot_labels(unit,"Management")
    b,_=stable_plot_labels(unit,"Grazing")
    output["robustness"]["joint_stable_management_and_grazing"]={
        "model":choose_fit(common,a&b,"joint_stable_management_and_grazing")
    }
    return output

def main():
    p=argparse.ArgumentParser()
    for n in ("npms","cover","receipt","eiv","out"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    ans=investigate(a.npms,a.cover,a.receipt,a.eiv)
    out=Path(a.out);out.parent.mkdir(exist_ok=True,parents=True)
    out.write_text(json.dumps(ans,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("NPMS_OBSERVATION_SOURCE_CLASSIFICATION",json.dumps({
        "source_units":ans["observed_survey_units"],
        "original_values":ans["source_aggregate_captions"],
        "source_coverage":ans["source_attribute_coverage"],
        "stability_results":{k:{
             "plots":v.get("stability",{}).get("plots_stable_all_years_labeled"),
             "status":v.get("model",v).get("status"),
             "beta":v.get("model",v).get("beta")
             } for k,v in ans["robustness"].items()},
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
