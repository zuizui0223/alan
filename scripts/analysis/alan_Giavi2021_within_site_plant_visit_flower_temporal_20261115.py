#!/usr/bin/env python3
"""Swiss ALAN true same-sampling-event floral cover x diurnal visitation timing.

Original Giavi et al 2021 Zenodo 4540407 Rdata MD5 f60b21d82c163cbb53520fda0bc61985.
Use site×plant seasonal sampling-event profiles on SAME observed flower opportunity.
Six paired independent site assignments, 37 sampled event pairs. The repeat field
is a survey index, not guaranteed equal-day spacing: overlap metric does NOT
need ordinal ordering. Distinct from pollen transfer or true visitor-flower
hourly synchrony / island shift. Do not infer unlisted plants were absent.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_Giavi2021_Swiss_original_daytime_guild_reanalysis_20261114 import (
    load_original_data, check_and_prepare, require
)

GUILDS=["All3","Coleoptera","Diptera","Hymenoptera"]
SUPPORT=[2,3,4]
SOURCE_DOI="10.1038/s41467-021-22011-8"
DATA_DOI="10.5281/zenodo.4540407"

def normalized_event_overlap(flower,visits):
    """Overlap of normalized nonnegative distributions on matched events.
    Only defined if flowers AND recorded visits >0, in BOTH light arms.
    """
    f=np.asarray(flower,dtype=float)
    v=np.asarray(visits,dtype=float)
    if len(f)!=len(v) or len(f)<2 or not np.isfinite(f).all() or not np.isfinite(v).all():
        raise ValueError("Mismatched/nonfinite/sparse common events")
    if (f<0).any() or (v<0).any():raise ValueError("Negative original counts")
    if f.sum()<=0 or v.sum()<=0:return None
    return float(np.minimum(f/f.sum(),v/v.sum()).sum())

def flower_proportional_null_expected_overlap(flower,visits,seed=20261115,B=2000):
    """Finite-visitation-count null: allocate total recorded visits across
    actual common events in proportion to the same arm's observed flower cover.
    This only subtracts a count/cover-dependent expected chance overlap:
    it does NOT randomize ALAN, test biological pollen mismatch, or restore
    unobserved plant × treatment rows.
    """
    f=np.asarray(flower,float); v=np.asarray(visits,float)
    if len(f)!=len(v) or len(f)<2 or not np.isfinite(f).all() or not np.isfinite(v).all():
        raise ValueError("Source event-count mismatch or nonfinite")
    if (f<0).any() or (v<0).any() or B<=0:
        raise ValueError("Bad original positive floral or visit counts")
    if f.sum()<=0 or v.sum()<=0:return None
    if not np.allclose(v,np.round(v),atol=1e-8):
        raise ValueError("Cannot multinomial-model source visit records with non-integer visit units")
    n=int(round(float(v.sum())));p=f/f.sum()
    draw=np.random.default_rng(seed).multinomial(n,p,size=B)
    return float(np.minimum(p,draw/n).sum(axis=1).mean())

def exact_site_sign_flip(site_deltas):
    """Exhaust all 2^n flips of true independently paired site summaries.
    Randomization schedule is unknown: diagnostic symmetry, NOT a formal
    randomized-assignment test and NOT multiple-test corrected.
    """
    import itertools
    v=np.asarray(site_deltas,float)
    if not len(v) or len(v)>12 or not np.isfinite(v).all():
        raise ValueError("Requires 1–12 finite site effects")
    obs=abs(v.mean())
    null=np.array([abs(np.mean(v*np.asarray(s)))
                   for s in itertools.product((-1.,1.),repeat=len(v))])
    return {
        "independent_site_count":int(len(v)),
        "two_sided_site_signflip_fraction_DESCRIPTIVE":
            float((null>=obs-1e-12).mean()),
        "NOT_formal_randomization_under_field_assigned_light":True,
        "nominal_sample_size_not_plant_site_time_rows":True,
    }

def holm_by_family(pvalues):
    """Monotone Holm adjustment across FOUR source guilds at k>=3 events."""
    order=sorted(pvalues,key=lambda k:pvalues[k])
    m=len(order);result={};last=0
    for i,k in enumerate(order):
        raw=float(pvalues[k])
        if not 0<=raw<=1:raise ValueError("P must be a valid probability")
        last=max(last,min(1.0,(m-i)*raw))
        result[k]=last
    return result

def site_bootstrap(delta,seed=20261115,iterations=10000):
    """Only six original experiment paired sites; never resample plant rows."""
    s=np.asarray(delta,float)
    if len(s)<3 or not np.isfinite(s).all():return None
    rng=np.random.default_rng(seed)
    sims=rng.choice(s,(iterations,len(s)),replace=True).mean(axis=1)
    return np.quantile(sims,[.025,.975]).tolist()

def source_grouping(a):
    a=check_and_prepare(a)
    # plab is duplicated once for each insect-order row. Preserve exact
    # plant event and treatment support before deduplicating abundance.
    keys=["site","pair","plant","treatment"]
    g=a.groupby(keys,observed=True).agg(
       floral_cover_min=("plab","min"),floral_cover_max=("plab","max"),
       floral_cover=("plab","median"),all3_visits=("dVis","sum"),
       guild_rows=("insect_group","nunique"),
       survey_repeat_min=("repeate","min"),survey_repeat_max=("repeate","max")
       ).reset_index()
    require(len(g)==628,"Published 628 plant x paired occasion x treatment cover observations")
    require(g.guild_rows.eq(3).all(),"Each observed plant event must have 3 measured insect orders")
    require(np.allclose(g.floral_cover_min,g.floral_cover_max),"3 duplicated source cover records differ")
    require(g.survey_repeat_min.eq(g.survey_repeat_max).all(),"Source repeat labels differ between orders")
    c=a.pivot_table(index=keys,columns="insect_group",values="dVis",
                    aggfunc="sum",observed=True).reset_index()
    require(len(c)==628,"Guild pivot unexpectedly lost source observed rows")
    g=g.merge(c,on=keys,validate="one_to_one")
    require(len(g)==628,"Source plant-event coverage not 628")
    require(g.pair.nunique()==37 and g.site.nunique()==6 and g.plant.nunique()==21,
            "Source site×plant×occasion scope changed")
    require(g.groupby("pair").site.nunique().eq(1).all(),
            "Source paired occasion spans multiple field sites")
    # Captures the number of true pairs at each experimental site, not
    # simply every positive flower record.
    design={
      "source_interaction_rows":int(len(a)),
      "original_source_plant_cover_event_rows":int(len(g)),
      "true_site_pairs":int(g.site.nunique()),
      "true_paired_sampling_occasions":int(g.pair.nunique()),
      "original_plant_species_labels":int(g.plant.nunique()),
      "source_cover_zero_records":int(g.floral_cover.eq(0).sum()),
      "source_positive_cover_records":int(g.floral_cover.gt(0).sum()),
      "source_zero_visits_at_observed_cover_records":int(g.all3_visits.eq(0).sum()),
      "source_diel_time_hour_not_available":True,
      "source_repeate_is_NOT_assumed_exact_calendar_day":True,
      "n_source_repeate_labels":int(g.survey_repeat_min.nunique()),
      "source_repeate_labels":sorted(str(v) for v in g.survey_repeat_min.unique()),
      "n_source_plant_observations_by_treatment":{
         str(k):int(v) for k,v in g.treatment.value_counts().items()},
      "paired_sampling_occasions_per_true_site":{
         "min":int(g.groupby("site").pair.nunique().min()),
         "max":int(g.groupby("site").pair.nunique().max())},
      "author_readme_required_for_repeate_calendar_mapping":True
    }
    return g,design

def observed_common_events(g):
    # Need flower cover actually measured in BOTH lit T and dark C at same
    # source site, paired sampling occasion and same plant species.
    idx=["site","pair","plant"]
    vals=["floral_cover","all3_visits","Coleoptera","Diptera","Hymenoptera"]
    wide=g.pivot(index=idx,columns="treatment",values=vals)
    # Keep only events with exact treatment opportunity values.
    both=wide.dropna(subset=[(f,arm) for f in vals for arm in ("C","T")])
    both.columns=[f"{v}_{a}" for v,a in both.columns]
    both=both.reset_index()
    require(len(both)<=314,"Too many same-event verified matched source plant rows")
    return both

def pair_plants(common,visit_col,support):
    """Pair within site×plant: only shared survey occasions, no missing=0."""
    source=[]
    for (site,plant),q in common.groupby(["site","plant"],observed=True):
        if len(q)<support:continue
        if q.pair.nunique()<support:continue
        t={}
        tnull={}
        for arm in ("C","T"):
            cov=q[f"floral_cover_{arm}"].to_numpy(float)
            visits=q[f"{visit_col}_{arm}"].to_numpy(float)
            score=normalized_event_overlap(cov,visits)
            if score is None:break
            # Deterministic, source-specific independent draw streams:
            # retained only as aggregate, never serialize plant/site keys.
            stream=int.from_bytes(__import__("hashlib").sha256(
                f"{site}/{plant}/{visit_col}/{support}/{arm}".encode()).digest()[:4],"big")
            null=flower_proportional_null_expected_overlap(cov,visits,seed=stream)
            t[arm]=score
            tnull[arm]=null
        if len(t)!=2:continue
        source.append({"source_site":site,"plant":plant,"source_common_occasion_n":len(q),
                       "lit_overlap":t["T"],"dark_overlap":t["C"],
                       "lit_minus_dark_temporal_overlap":t["T"]-t["C"],
                       "lit_minus_dark_null_expected_temporal_overlap":tnull["T"]-tnull["C"],
                       "lit_minus_dark_excess_overlap_beyond_visit_count_null":
                           (t["T"]-tnull["T"])-(t["C"]-tnull["C"]),
                       "lit_total_visits_same_support":float(q[f"{visit_col}_T"].sum()),
                       "dark_total_visits_same_support":float(q[f"{visit_col}_C"].sum())})
    return source

def summarize_paired(source,seed):
    if not source:return {"supported_plant_site_series":0,"supported_independent_sites":0}
    z=pd.DataFrame(source)
    bysite=z.groupby("source_site").lit_minus_dark_temporal_overlap.median().sort_index()
    d=bysite.to_numpy(float)
    adj=z.groupby("source_site").lit_minus_dark_excess_overlap_beyond_visit_count_null.median().sort_index()
    d_adj=adj.to_numpy(float)
    sample_bias=z.groupby("source_site").lit_minus_dark_null_expected_temporal_overlap.median().sort_index().to_numpy(float)
    # Output no source site or rare plant identity.
    return {
      "supported_plant_site_series":len(z),
      "supported_plant_taxa":int(z.plant.nunique()),
      "supported_independent_sites":len(bysite),
      "true_min_source_occasions_per_series":int(z.source_common_occasion_n.min()),
      "true_max_source_occasions_per_series":int(z.source_common_occasion_n.max()),
      "observed_plant_site_positive_delta_n":int((z.lit_minus_dark_temporal_overlap>0).sum()),
      "observed_plant_site_negative_delta_n":int((z.lit_minus_dark_temporal_overlap<0).sum()),
      "mean_of_site_medians_lit_minus_dark_overlap":float(d.mean()),
      "median_of_site_medians_lit_minus_dark_overlap":float(np.median(d)),
      "site_bootstrap95_mean_lit_minus_dark_DESCRIPTIVE":site_bootstrap(d,seed),
      "source_sites_positive_site_median_n":int((d>0).sum()),
      "source_sites_negative_site_median_n":int((d<0).sum()),
      "source_site_level_effects_anonymized":[float(v) for v in d],
      "observed_raw_overlap_site_signflip_diagnostic":exact_site_sign_flip(d),
      "observed_site_equal_overlaps":{
          "lit_site_median_mean":float(z.groupby("source_site").lit_overlap.median().mean()),
          "dark_site_median_mean":float(z.groupby("source_site").dark_overlap.median().mean())
      },
      "visit_count_and_cover_proportional_null_expected_delta_site_equal":float(sample_bias.mean()),
      "lit_minus_dark_null_excess_event_overlap_site_equal":float(d_adj.mean()),
      "lit_minus_dark_null_excess_site_bootstrap95_DESCRIPTIVE":site_bootstrap(d_adj,seed+40000),
      "site_median_null_excess_positive_n":int((d_adj>0).sum()),
      "site_median_null_excess_negative_n":int((d_adj<0).sum()),
      "null_excess_site_level_anonymized":[float(v) for v in d_adj],
      "null_excess_site_signflip_diagnostic":exact_site_sign_flip(d_adj),
      "null_assumption":"Conditional on original visit total and observed same-arm flower cover across shared events, assign visits by independent multinomial proportional to cover. 2000 draws per original site-plant-arm. NOT a randomized lighting test.",
      "null_mc_draws_each_arm":2000,
      "no_effective_pollen_or_original_causal_randomization_recovered":True,
      "coverage_conditional_on_visitor_records_in_both_arms":True,
      "source_missing_excluded_not_zero":True,
      "not_a_pollinator_efficacy_or_true_phenological_mismatch":True
    }

def run():
    a,receipt=load_original_data()
    g,design=source_grouping(a)
    common=observed_common_events(g)
    design.update({
      "plant_observations_both_treatments_on_same_site_pair":int(len(common)),
      "source_unique_plant_site_common_occasion_groups":
          int(common.groupby(["site","plant"]).ngroups),
      "source_pair_x_plant_unmatched_rows":int(len(g)-2*len(common)),
      "n_common_occasions_per_plant_site":{
          "min":int(common.groupby(["site","plant"]).size().min()),
          "max":int(common.groupby(["site","plant"]).size().max())},
      "unmatched_source_rows_are_not_coded_biological_absence":True
    })
    compare={}
    for i,visit in enumerate(GUILDS):
        name="all_three_documented_diurnal_orders" if visit=="All3" else visit
        c="all3_visits" if visit=="All3" else visit
        compare[name]={}
        for k in SUPPORT:
            series=pair_plants(common,c,k)
            compare[name][str(k)]=summarize_paired(series,20261115+i*100+k)
    # Four source groups are an explicit exploratory multiplicity family.
    primary={k:v["3"].get("null_excess_site_signflip_diagnostic",{})
             .get("two_sided_site_signflip_fraction_DESCRIPTIVE")
             for k,v in compare.items()
             if v["3"].get("supported_independent_sites",0)>=3}
    family_adjusted=holm_by_family(primary)
    return {
      "four_guild_primary_threshold3_exact_site_signflip_Holm_exploratory":family_adjusted,
      "study_family_k":4,
      "multiple_thresholds_2_3_4_not_formally_adjusted":True,
      "schema":"ALAN_Swiss_Giavi2021_original_same_plant_event_diurnal_flower_visitor_overlap_v1",
      "study_DOI":SOURCE_DOI,"original_data_DOI":DATA_DOI,
      "source_Rdata_md5":receipt["original_md5"],
      "source_Rdata_sha256":receipt["sha256"],
      "original_source_readme_author_design_not_fake_sample_date":True,
      "original_data_support":design,
      "plant_flower_cover_vs_visits_matched_temporal_overlap":compare,
      "main_result_exploratory_study_level_only":True,
      "no_hourly_open_flower_vs_nocturnal_moth_profile_available":True,
      "no_fruit_seed_or_effective_pollen_measured":True,
      "no_geographic_island_contrast_measured":True,
      "missing_unobserved_flower_source_records_not_imputed_zero":True,
      "no_crosspaper_Grognuz2026_site_date_join":True,
      "interpretation":"Temporal congruence across repeated field survey EVENTS for co-observed flowers and diurnal visitors; not circadian flower-pollinator mismatch and not causal fitness."
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    args=p.parse_args();out=Path(args.out)
    if out.exists():raise FileExistsError("Do not overwrite")
    z=run()
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(z,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_GIAVI2021_TRUE_FLORAL_VISITOR_TEMPORAL_SUPPORT",json.dumps({
      "source":z["original_data_support"],
      "results":z["plant_flower_cover_vs_visits_matched_temporal_overlap"]
      },ensure_ascii=False))
if __name__=="__main__":main()
