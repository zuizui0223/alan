#!/usr/bin/env python3
"""Robustness of 39-island ALAN × externally measured diel reproductive function.

Source-locked full93m PHOTO IDENTITIES, not population counts or observed
pollen deposition. Two direct checks:
  A) Leave one photographed experimental species / island out of the fixed,
     original within-island-year co-documented night-vs-day observed-flora
     regression. No cherry-picked omission or altered inferential claim.
  B) Independently verified species-direct autonomous selfing and SC/SI
     trait overlap with externally studied night-vs-day reproductive effects.

Missing published evidence is UNKNOWN, never a biological zero.
All outputs are aggregate-only, not restricted plant traits or coordinates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from scripts.analysis.alan_real_daynight_reproductive_PCNL_562_20261113 import (
    experimental_daynight, KEY, LIGHT, coobserved_functional_ratio,
)
from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    load_sites, load_labels, original_accepted_species, site_outcome_model, require,
)

LEDGER_SHA="7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577"
BASELINE=-1.1747574942083159
N_ISLANDS=562
N_SITES=77647
N_SPECIES_SITE=2012925


def original_counts(d, exposure, sites):
    relevant=exposure.loc[
        exposure.night_stronger|exposure.day_stronger,
        ["species", "night_stronger", "day_stronger", "night_share","n_studies"]
    ].copy()
    source=d.loc[d.species.isin(relevant.species),KEY+["species"]].copy()
    source=source.merge(relevant,on="species",how="left",validate="many_to_one")
    require(source[["night_stronger","day_stronger"]].notna().all().all(),
            "Original photo species failure to join known author experimental effects")
    source["genus"]=source.species.str.split().str[0]
    aggregate=(source.groupby(KEY,observed=True).agg(
        night=("night_stronger","sum"),day=("day_stronger","sum"))
        .reset_index())
    s=sites.merge(aggregate,on=KEY,how="left",validate="one_to_one")
    for label in ("night","day"):s[label]=s[label].fillna(0).astype(int)
    require((s.night+s.day<=s.all_accepted_recorded_species).all(),
            "Experimental source taxon pair double counted")
    require(int((s.night>=1).sum())>0 and int((s.day>=1).sum())>0,
            "No original documented functional channels in photo sources")
    require(int(((s.night>=1)&(s.day>=1)).sum()) == 2565,
            "Original 2565 co-documented original photo-cell support changed")
    # Unique taxon x source photographic cell/year row for omission audit.
    require(not source.duplicated(KEY+["species"]).any(), "Repeated source plant")
    return s,source,relevant


def supported(s):
    """Require BOTH reproductive groups at each photo cell, then spatially
    simultaneous two illuminated sites with varying previous light per
    actual island x year (same gate as original PR46)"""
    eligible=s.loc[s.night.ge(1)&s.day.ge(1)].copy()
    g=eligible.groupby(["island","year"],observed=True)[LIGHT].agg(["size","min","max"])
    valid=g.loc[g["size"].ge(2)&(g["max"]-g["min"]).ge(.02)]
    q=eligible.merge(valid[[]].reset_index(),on=["island","year"],
                     how="inner",validate="many_to_one")
    return q,dict(positive_sites_both=int(len(eligible)),
                  positive_islands_both=int(eligible.island.nunique()),
                  eligible_island_years=int(len(valid)),
                  evaluated_sites=int(len(q)),
                  evaluated_islands=int(q.island.nunique()))


def fit_sample(s, *, min_islands=25):
    q,support=supported(s)
    if len(q)<150 or q.island.nunique()<min_islands:
        return {"status":"not_estimable_common_source_geography","support":support}
    ratio=q.night.to_numpy(float)/(q.night+q.day).to_numpy(float)
    raw=site_outcome_model(q,ratio*q.all_accepted_recorded_species.to_numpy(float))
    return {"status":"source_photographic_composition_observational",
            "support":support, "estimate":raw}


def get_beta(result):
    return result.get("estimate",{}).get("beta_pp")


def influence(s,source):
    baseline=fit_sample(s)
    require(baseline["status"]=="source_photographic_composition_observational",
            "Original source coobservation not estimable")
    require(baseline["support"]["evaluated_sites"]==2373 and
            baseline["support"]["evaluated_islands"]==39,
            "Original 39 islands/2373 photo cells not recovered")
    require(abs(get_beta(baseline)-BASELINE)<1e-8,
            "Original source coefficient drift")

    original_sites,_=supported(s)
    islands=sorted(original_sites.island.unique())
    leave_island=[]
    for island in islands:
        z=original_sites.loc[original_sites.island.ne(island)].copy()
        # Every retained island has the same source-supported island×year
        # sites and PCNL range. No reselecting days/time/spatial thresholds.
        ratio=z.night.to_numpy(float)/(z.night+z.day).to_numpy(float)
        x=site_outcome_model(z,ratio*z.all_accepted_recorded_species.to_numpy(float))
        leave_island.append({"left_out_island":str(island),
                             "beta":x["beta_pp"],
                             "ci":x["95CI_pp"],
                             "p_nominal":x["nominal_p"]})
    # Public aggregate must NOT expose island IDs for rare study taxa.
    island_betas=np.array([z["beta"] for z in leave_island])
    island_ci=[z["ci"] for z in leave_island]

    # Source is only ~20,000 positive study-specific species×site-year rows.
    # Remove one experimental PHOTO TAXON, not a study's raw observation row.
    site_index=pd.MultiIndex.from_frame(s[KEY])
    positions=site_index.get_indexer(pd.MultiIndex.from_frame(source[KEY]))
    require((positions>=0).all(),"Original source plant key not in 562 island frame")
    source=source.assign(_site_ix=positions)
    omission=[]
    for name,part in source.groupby("species",observed=True):
        group=part[["night_stronger","day_stronger"]].iloc[0]
        require(group.night_stronger != group.day_stronger,
                "Species treatment class ambiguous")
        changed=s.copy()
        col="night" if bool(group.night_stronger) else "day"
        src=changed[col].to_numpy(copy=True)
        src[part._site_ix.to_numpy(int)]-=1
        require((src>=0).all(),"Source plant present more than once at a site")
        changed[col]=src
        est=fit_sample(changed)
        omission.append({
            "source_photo_site_years_for_taxon":len(part),
            "study_count":int(part.n_studies.iloc[0]),
            "species_functional_class":col,
            "status":est["status"],
            "source_islands_after_omission":int(est["support"].get("evaluated_islands",0)),
            "source_sites_after_omission":int(est["support"].get("evaluated_sites",0)),
            "beta":get_beta(est),
            "ci":est.get("estimate",{}).get("95CI_pp"),
            "p_nominal":est.get("estimate",{}).get("nominal_p"),
            # NO species name or raw species x coordinate exported.
        })
    supported_omissions=[z for z in omission if z["beta"] is not None]
    values=np.array([z["beta"] for z in supported_omissions])
    def stats(vals):
        vals=np.asarray(vals,float)
        return {"n":len(vals), "range":[float(vals.min()),float(vals.max())],
                "median":float(np.median(vals)),
                "fraction_negative":float((vals<0).mean()),
                "n_CI_entirely_below_zero":int(sum(
                    z["ci"][1]<0 for z in supported_omissions))}

    summary={
        "baseline_39island_codocumented":baseline,
        "leave_one_island_out":{
            "island_deletions":len(leave_island),
            "beta_range":[float(island_betas.min()),float(island_betas.max())],
            "median_beta":float(np.median(island_betas)),
            "fraction_negative":float((island_betas<0).mean()),
            "n_nominal_95CI_excludes_zero":int(sum(c[1]<0 or c[0]>0 for c in island_ci)),
            "n_beta_sign_reversed":int(sum(v>=0 for v in island_betas)),
            "all_island_identities_suppressed":True,
        },
        "leave_one_published_experiment_species_out":{
            "total_experimental_species_deletions":len(omission),
            "supported_after_deletion":len(supported_omissions),
            "not_estimable_after_deletion":len(omission)-len(supported_omissions),
            "beta_distribution":stats(values) if len(values) else None,
            "n_source_site_below_original":int(sum(
                z["source_sites_after_omission"]<baseline["support"]["evaluated_sites"]
                for z in supported_omissions)),
            "largest_abs_coefficient_delta_top_5_without_identity":[
                {k:z[k] for k in ("species_functional_class",
                       "source_photo_site_years_for_taxon",
                       "study_count","source_islands_after_omission",
                       "source_sites_after_omission","beta")}
                for z in sorted(supported_omissions,
                    key=lambda x:abs(x["beta"]-BASELINE),reverse=True)[:5]
            ],
            "no_taxon_names_and_no_sites_exported":True,
        }
    }
    return summary



def time_and_genus_sensitivities(s,source,exposure):
    """Predeclared robustness: whole-genus deletion, source time normalization,
    measured night-share continuous on SAME 39-island cohort, site weighting."""
    from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import fixed_effects_regression
    index=pd.MultiIndex.from_frame(s[KEY])
    src=source.copy()
    src["_site_ix"]=index.get_indexer(pd.MultiIndex.from_frame(src[KEY]))
    require(src._site_ix.ge(0).all(),"Source species outside original photo cells")
    deletion=[]
    for genus,rows in src.groupby("genus",observed=True):
        v=s.copy()
        for label,flag in [("night","night_stronger"),("day","day_stronger")]:
            keep=rows.loc[rows[flag]]
            arr=v[label].to_numpy(copy=True)
            arr[keep._site_ix.to_numpy(int)]-=1
            require((arr>=0).all(),"Genus omission double subtract")
            v[label]=arr
        result=fit_sample(v)
        deletion.append({
            "genus_source_taxa":int(rows.species.nunique()),
            "source_species_site_year_rows":len(rows),
            "status":result["status"],
            "beta":get_beta(result),
            "CI95":result.get("estimate",{}).get("95CI_pp"),
            "remaining_islands":result["support"].get("evaluated_islands"),
            "remaining_source_sites":result["support"].get("evaluated_sites"),
        })
    valid=[d for d in deletion if d["beta"] is not None]
    require(valid,"No supported genus deletion for source ratio")
    vals=np.array([d["beta"] for d in valid])
    q,base=supported(s)
    require(len(q)==2373 and q.island.nunique()==39,"Baseline sites changed")

    # Same 39 islands/site-years; estimate continuous original external study
    # night share, conditional on BOTH experimentally classified plant classes.
    grouped=src.groupby(KEY,observed=True).agg(
        n_trait_species=("species","size"),
        sum_night_share=("night_share","sum")).reset_index()
    con=q.merge(grouped,on=KEY,validate="one_to_one")
    require(con.n_trait_species.ge(2).all(),"Each source co-site must have both classes")
    continuous=con.sum_night_share.to_numpy(float)/con.n_trait_species.to_numpy(float)
    fit_cont=site_outcome_model(con,continuous*con.all_accepted_recorded_species.to_numpy(float))

    # Same light/same-site groups, only weight changes: decompose whether
    # populous documented islands dominate the sign.
    ix=q.island.astype(str).str.cat(q.year.astype(str),sep="|").to_numpy()
    fit_site_equal=fixed_effects_regression(
        (q.night.to_numpy(float)/(q.night+q.day).to_numpy(float)),
        q[LIGHT].to_numpy(float),q.log_obs.to_numpy(float),
        np.ones(len(q)),q.island.to_numpy(),[ix])

    # Alternate definition: normalize laboratory/field daytime and nighttime
    # exclusion-response by eligible pollination-hours, no mean-time imputation.
    hours=src.merge(exposure[["species","night_share_per_hour"]],
                    on="species",validate="many_to_one")
    hours["hour_night"]=hours.night_share_per_hour.gt(.5+1e-12)
    hours["hour_day"]=hours.night_share_per_hour.lt(.5-1e-12)
    agg=hours.groupby(KEY,observed=True)[["hour_night","hour_day"]].sum().reset_index()
    alt=s[KEY+list(set(s.columns)-set(KEY))].merge(agg,on=KEY,how="left",
                                         validate="one_to_one")
    alt["night"]=alt.hour_night.fillna(0).astype(int)
    alt["day"]=alt.hour_day.fillna(0).astype(int)
    altfit=fit_sample(alt)
    return {
        "leave_one_genus_out":{
            "total_studied_genera":len(deletion),
            "supported_genus_deletions":len(valid),
            "not_estimable_genus_deletions":len(deletion)-len(valid),
            "beta_range":[float(vals.min()),float(vals.max())],
            "median_beta":float(np.median(vals)),
            "fraction_negative":float((vals<0).mean()),
            "n_nominal_95CI_entirely_negative":int(sum(z["CI95"][1]<0 for z in valid)),
            "n_beta_sign_reversed":int(sum(z["beta"]>=0 for z in valid)),
            "largest_genus_omission_effect_no_identity":sorted(
                [{k:z[k] for k in ("genus_source_taxa","source_species_site_year_rows",
                     "beta","remaining_islands","remaining_source_sites")}
                 for z in valid],key=lambda a:abs(a["beta"]-BASELINE),reverse=True)[:5],
            "no_taxonomic_or_site_identity_released":True,
        },
        "same_39island_original_coobservation_cohort":{
            "original_source_photo_cell_years":len(q),
            "original_islands":int(q.island.nunique()),
            "continuous_night_share_from_external_experiment":fit_cont,
            "site_year_equal_instead_of_island_equal_bin_class":fit_site_equal,
        },
        "hours_normalized_alternate_functional_classification":{
            "observed_reproductive_photo_taxa_with_valid_hours":int(
                exposure.loc[exposure.species.isin(src.species)
                             & exposure.night_share_per_hour.notna(),"species"].nunique()),
            "alternate_class_source_support":altfit,
            "source_group_is_altered_by_day_night_hours_missingness":True,
        },
        "interpretation_limit":"None of these are same-site experimental pollen transfer, compensating pollination, or an independent replication",
    }


def direct_selfing(exposure,source,ledger:Path):
    require(hashlib.sha256(ledger.read_bytes()).hexdigest()==LEDGER_SHA,
            "Original island species direct trait file changed")
    traits=pd.read_csv(ledger, dtype=str,keep_default_na=False)
    require((len(traits),traits.accepted_species.nunique())==(134268,62423),
            "Original published direct botanical traits changed")
    candidates=traits.loc[traits.trait_name.isin([
        "autonomous_selfing_capacity","self_incompatibility","mating_system"])
        & traits.resolution_status.eq("resolved")
        & traits.quality.isin(["high","medium"])
        & traits.multistate.eq("False")
        & traits.normalized_value.ne("")
        & ~traits.normalized_value.str.contains("|",regex=False)].copy()
    require(not candidates.duplicated(["accepted_species","trait_name"]).any(),
            "Trait source conflicts")
    t=candidates.pivot(index="accepted_species",columns="trait_name",
                       values="normalized_value")
    if t.empty:
        return {"status":"not_estimable_no_direct_selfing_values"}
    study=exposure.copy()
    study=study.merge(t,left_on="species",right_index=True,how="left",
                      validate="one_to_one")
    photo_names=set(source.species)
    study["in_original_photo_cohort"]=study.species.isin(photo_names)
    study["genus"]=study.species.str.split().str[0]

    result={"source_original_direct_trait_taxa":int(t.index.nunique()),
            "source_exp_reproductive_species":len(study),
            "source_photo_matched_experiment_species":int(study.in_original_photo_cohort.sum()),
            "tiers":{}}
    for scope in ("all_134_experimental_species","photo_matched_52_experimental_species"):
        d=study if scope.startswith("all_") else study[study.in_original_photo_cohort]
        present=d.night_stronger; absent=d.day_stronger
        require((present|absent).all(),"Unknown experiment binary direction unexpectedly included")
        kinds={}
        for trait in ("autonomous_selfing_capacity","self_incompatibility","mating_system"):
            vals=d[trait] if trait in d else pd.Series(pd.NA,index=d.index)
            known=vals.notna()
            valid=vals.isin(["autonomous","absent"] if trait=="autonomous_selfing_capacity" else
                            ["SC","SI"] if trait=="self_incompatibility" else
                            ["predominantly_selfing","obligate_selfing","predominantly_outcrossing","mixed_mating"])
            q=d.loc[valid]
            stats={
                "trait_known_any":int(known.sum()),
                "trait_known_valid":int(valid.sum()),
                "night_higher_trait_known":int((present&valid).sum()),
                "day_higher_trait_known":int((absent&valid).sum()),
                "quality_source":"direct_species_high_medium",
                "source_trait_category_counts":vals.value_counts().to_dict()
            }
            # Prespecified 2x2 only if > 15 independent taxa, and >=4
            # plants in each observed trait category. Otherwise no p value.
            if trait=="autonomous_selfing_capacity":
                positive="autonomous"; negative="absent"
            elif trait=="self_incompatibility":
                positive="SC"; negative="SI"
            else:
                positive="predominantly_selfing"; negative="predominantly_outcrossing"
            nplus=int(q[trait].eq(positive).sum()); nminus=int(q[trait].eq(negative).sum())
            a=int((q.night_stronger&q[trait].eq(positive)).sum())
            b=int((q.night_stronger&q[trait].eq(negative)).sum())
            if nplus>=4 and nminus>=4 and (nplus+nminus)>=15:
                odds,p=fisher_exact([[a,b],[nplus-a,nminus-b]])
                stats["source_trait_night_association"]=dict(
                    trait_positive=positive,trait_negative=negative,
                    positive_total=nplus,negative_total=nminus,
                    night_positive=a,night_negative=b,
                    odds_ratio=float(odds) if np.isfinite(odds) else None,
                    exploratory_p=float(p))
            else:
                stats["source_trait_night_association"]={
                    "status":"NOT_ESTIMABLE_INDEPENDENT_SPECIES_SUPPORT",
                    "positive_total":nplus,"negative_total":nminus}
            kinds[trait]=stats
        result["tiers"][scope]=kinds
    result["inferential_ceiling"]=[
        "Self-compatibility is not autonomous selfing",
        "Species-direct trait observed externally not confirmed at photographed island",
        "Nighttime treatment advantage is not moth identity or effective pollen service",
        "No experimental proof of ecological resilience or causal light impact",
        "Do not infer non-selfing from missing source"
    ]
    return result


def run(sites_path,globisource,groups,study,traitfile,ledger,out):
    out=Path(out)
    if out.exists():raise FileExistsError(out)
    sites=load_sites(sites_path)
    require((len(sites),sites.island.nunique())==(N_SITES,N_ISLANDS),
            "Original 562-island frame changed")
    old=load_labels(globisource)
    flora=original_accepted_species(groups,sites,old)
    require(len(flora)==N_SPECIES_SITE,"Original species×site source changed")
    experiment, meta=experimental_daynight(study,traitfile)
    siteframe,source,_=original_counts(flora,experiment,sites)
    source_stability=influence(siteframe,source)
    selfing=direct_selfing(experiment,source,ledger)
    third_axis=time_and_genus_sensitivities(siteframe,source,experiment)
    output={
       "schema":"ALAN_39island_codocumented_stability_selfing_direct_v1",
       "status":"AUTHENTIC_FROZEN_SOURCE_REPLAY_EXPLORATORY",
       "original_cohort":{
         "global_islands":int(sites.island.nunique()),
         "source_site_years":len(sites),
         "source_accepted_species":int(flora.species.nunique()),
         "source_plant_site_years":len(flora),
         "PCNL":"1992 to prior year corrected annual light, mean log1p",
         "recorded_time":"original 2019/2021/2023/2024/2025 plant observation years"},
       "external_experiment_source":meta,
       "functional_coobservation_influence_diagnostics":source_stability,
       "external_autonomous_selfing_and_SCSI_species_direct_overlap":selfing,
       "site_weight_and_genus_time_normalized_sensitivity":third_axis,
       "guardrails":{
         "no_anonymized_site_x_species_location_data_exported":True,
         "unknown_source_function_not_absence":True,
         "study_selection_not_biological_planted_census":True,
         "source_replication_not_causal_ALAN_effect":True,
         "classification_inconsistent_across_studies_not_confirmatory":True,
         "study_selected_p_values_posthoc":True}
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_39ISLAND_LEAVEOUT_AND_DIRECT_SELFING",
          json.dumps({"stability":source_stability,
                      "selfing":selfing},ensure_ascii=False),flush=True)
    return output


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for n in ["sites","globi","groups","study","traitfile","ledger","out"]:
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.study,a.traitfile,a.ledger,a.out)
