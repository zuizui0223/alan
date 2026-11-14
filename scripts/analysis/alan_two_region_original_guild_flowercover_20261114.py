#!/usr/bin/env python3
"""Real Switzerland–Tibet independent field study contrast: guild counts vs floral opportunity.

Rehydrate both *genuine original* published field datasets, SHA/MD5-verify
authors' releases through the prior scripts, and compare 3 actual shared
daytime insect orders under light within their own study-specific units.

Swiss 2016 Giavi: six paired original meadow sites, 37 paired sampling
events, plab = measured flowering-plant cover proxy. Deduplicate plab
once per site × paired event × plant × treatment (NOT one copy per order).
Tibet 2024 Lu: five independently assigned original field blocks, four
sampling rounds, raw visitation counts and survey observation opportunities;
flower cover per plant is NOT available for direct normalization.

Report:
1. independent site/block-equal log visit ratios and full-site bootstrap CIs,
2. Swiss original flower-cover-standardized visit rate sensitivity,
3. 3-order Jensen-Shannon composition turnover per independent field unit,
4. fly versus bee relative-visitor *channel balance*, never pollen compensation,
5. within-Tibet night-vs-day contrasts and heterogeneity checks.

No claims of equal lighting spectra, global meta-analytic effectiveness,
new plant pollen deposition, phenological mismatch, predation, island ALAN
vulnerability or taxonomic guild absence from missing source records.
"""
from __future__ import annotations

import argparse, hashlib, io, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from openpyxl import load_workbook

from scripts.analysis.alan_Giavi2021_Swiss_original_daytime_guild_reanalysis_20261114 import (
    load_original_data as swiss_source, check_and_prepare as swiss_validate,
)
from scripts.analysis.alan_original_Lu2026_daynight_guild_syndrome_source_20261114 import (
    download as tibet_source, sheet as tibet_sheet,
    original_visitation_design as tibet_validate, XLSX_SHA,
)

ORDERS=("Diptera","Hymenoptera","Coleoptera")
LU={"Diptera":"D_abundance","Hymenoptera":"H_abundance",
    "Coleoptera":"C_abundance","Lepidoptera":"L_abundance",
    "Other":"O_abundance"}
PSEUDO=.5
SEED=20261114
NBOOT=5000

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_2REGION_REAL_COMMON_GUILDS: "+str(msg))

def lograte(n_lit,n_dark,expo_lit,expo_dark,pseudo=PSEUDO):
    require(min(n_lit,n_dark)>=0 and min(expo_lit,expo_dark)>0,"Unusable visit opportunity")
    return float(math.log((n_lit+pseudo)/(n_dark+pseudo)) -
                 math.log(expo_lit/expo_dark))

def jsbits(a,b):
    """Jensen–Shannon divergence (bits), valid even if an original guild has zero visits."""
    p=np.asarray(a,dtype=float);q=np.asarray(b,dtype=float)
    require(np.all(p>=0) and np.all(q>=0),"Negative original guild counts")
    require(p.sum()>0 and q.sum()>0,"No biological visitor community observed")
    p=p/p.sum();q=q/q.sum();m=(p+q)/2
    s=0.
    for x in (p,q):
        mask=x>0
        s+=.5*float(np.sum(x[mask]*np.log2(x[mask]/m[mask])))
    require(-1e-10<=s<=1+1e-10,"Impossible source guild divergence")
    return float(max(s,0))

def standardize_swiss():
    raw,receipt=swiss_source()
    d=swiss_validate(raw)
    require(len(d)==1884 and d.site.nunique()==6 and d.pair.nunique()==37,
            "Published Swiss 2016 original source changed")
    groups=d.groupby(["site","treatment","insect_group"],observed=True).dVis.sum()
    # Repeated plab is not three independent flower-cover opportunities.
    plab=d.groupby(["site","pair","plant","treatment"],observed=True).plab.agg(
         ["median","min","max"]).reset_index()
    incons=int((plab["max"]-plab["min"]>1e-9).sum())
    cover=plab.groupby(["site","treatment"]).median.sum().unstack("treatment")
    require(set(cover.columns)=={"C","T"} and cover.notna().all().all() and
            cover.gt(0).all().all(),"Original Swiss floral cover not observed in both treatments")
    site_ids=sorted(d.site.unique())
    records=[]
    for i,site in enumerate(site_ids,1):
        sitecounts={}
        for guild in ORDERS:
            c=float(groups.loc[(site,"C",guild)])
            t=float(groups.loc[(site,"T",guild)])
            sitecounts[guild]=(c,t)
        obs_c=float(d.loc[d.site.eq(site)&d.treatment.eq("C")].pair.nunique())
        obs_t=float(d.loc[d.site.eq(site)&d.treatment.eq("T")].pair.nunique())
        require(obs_c==obs_t and obs_c>=2,"Swiss field sampling events unbalanced")
        cover_c=float(cover.loc[site,"C"]);cover_t=float(cover.loc[site,"T"])
        row={"anonymous_original_site_index":i,
             "paired_survey_events":int(obs_c),
             "cover_dark":cover_c,"cover_lit":cover_t,
             "flower_cover_ratio_lit_over_dark":cover_t/cover_c,
             "guilds":{}}
        for guild,(c,t) in sitecounts.items():
            row["guilds"][guild]={
                "dark_original_visits":c,"lit_original_visits":t,
                "log_lit_dark_per_survey":lograte(t,c,obs_t,obs_c),
                "log_lit_dark_per_recorded_flower_cover":lograte(t,c,cover_t,cover_c)}
        row["three_order_JS_divergence_bits"]=jsbits(
            [sitecounts[k][0] for k in ORDERS],[sitecounts[k][1] for k in ORDERS])
        records.append(row)
    return records,{
        "original_data_Rdata_SHA256":receipt["sha256"],
        "field_experiments":6,"paired_survey_events":37,
        "flower_cover_records_deduplicated":int(len(plab)),
        "flower_cover_inconsistent_within_repeated_guilds":incons,
        "original_plant_taxa":21,"original_record_rows":len(d),
        "source_guild_visits":int(d.dVis.sum()),
        "Swiss_plab_is_discrete_flower_cover_proxy_not_flower_count":True,
        "Swiss_flower_cover_deduplicated_once_per_real_plant_event":True
    }

def standardize_tibet():
    raw,outer=tibet_source()
    w=load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    vis=tibet_validate(tibet_sheet(w,"visitation"))
    w.close()
    require(len(vis)==400 and vis.block.nunique()==5,"Original Lu five-block count drift")
    records={"Diurnal":[],"Nocturnal":[]}
    for period in records:
        d=vis.loc[vis.sampling_time.eq(period)].copy()
        for block in range(1,6):
            a=d.loc[d.block.eq(block)]
            n_dark=int(a.treatment.eq("Unlit").sum())
            n_lit=int(a.treatment.eq("Lit").sum())
            require(n_dark>=10 and n_lit>=10,
                    "Original Tibet survey block insufficient light/dark opportunities")
            c={}
            for k in ORDERS:
                c[k]=(float(a.loc[a.treatment.eq("Unlit"),LU[k]].sum()),
                      float(a.loc[a.treatment.eq("Lit"),LU[k]].sum()))
            row={"anonymous_original_field_block":block,
                  "survey_observation_opportunities_dark":n_dark,
                  "survey_observation_opportunities_lit":n_lit,"guilds":{}}
            for guild,(dark,lit) in c.items():
                row["guilds"][guild]={
                    "dark_original_visits":dark,"lit_original_visits":lit,
                    "log_lit_dark_per_survey":lograte(lit,dark,n_lit,n_dark)}
            row["three_order_JS_divergence_bits"]=jsbits(
                [c[k][0] for k in ORDERS],[c[k][1] for k in ORDERS])
            records[period].append(row)
    return records,{
        "source_original_workbook_SHA256":XLSX_SHA,
        "source_EuropePMC_wrapper_SHA_is_dynamic":True,
        "original_record_rows":len(vis),
        "field_experimental_blocks":5,
        "original_survey_rounds":4,
        "original_independent_regions":1,
        "true_flower_cover_available":False,
        "source_row_conservation_not_assumed_across_5_orders":True,
        "n_records_dark_day":int(vis.loc[
            vis.sampling_time.eq("Diurnal")&vis.treatment.eq("Unlit")].shape[0]),
        "n_records_lit_day":int(vis.loc[
            vis.sampling_time.eq("Diurnal")&vis.treatment.eq("Lit")].shape[0])
    }

def summarize(records,cover=False):
    col="log_lit_dark_per_recorded_flower_cover" if cover else "log_lit_dark_per_survey"
    rng=np.random.default_rng(SEED)
    out={}
    for guild in ORDERS:
        v=np.array([r["guilds"][guild][col] for r in records],dtype=float)
        require(len(v)>=5 and np.isfinite(v).all(),"Original independent site sample invalid")
        boot=rng.choice(v,size=(NBOOT,len(v)),replace=True).mean(axis=1)
        raw_dark=sum(r["guilds"][guild]["dark_original_visits"] for r in records)
        raw_lit=sum(r["guilds"][guild]["lit_original_visits"] for r in records)
        out[guild]={
            "independent_sites_or_blocks":len(v),
            "raw_original_dark_visits":raw_dark,
            "raw_original_lit_visits":raw_lit,
            "raw_pooled_lit_over_dark_ratio":float(raw_lit/raw_dark) if raw_dark else None,
            "equal_site_geometric_mean_lit_dark_VISIT_RATE_ratio":float(np.exp(v.mean())),
            "equal_site_log_visit_ratio":float(v.mean()),
            "descriptive_95CI_site_bootstrap_ratio":list(map(float,
                np.exp(np.quantile(boot,[.025,.975])))),
            "original_sites_where_lit_visit_rate_higher":int(np.sum(v>0)),
            "original_sites_where_lit_visit_rate_lower":int(np.sum(v<0)),
            "source_anonymous_site_ratios":[float(np.exp(x)) for x in v],
            "no_claimed_causal_pollen_service_or_global_meta_effect":True}
    div=np.array([r["three_order_JS_divergence_bits"] for r in records])
    out["three_order_JS_composition"]= {
        "source_n_independent_field_units":len(div),
        "mean_JS_divergence_bits":float(div.mean()),
        "median_JS_divergence_bits":float(np.median(div)),
        "source_5or6_unit_range":list(map(float,[div.min(),div.max()])),
        "visitation_composition_not_pollen_transport":True}
    v=np.array([r["guilds"]["Diptera"][col]-
                r["guilds"]["Hymenoptera"][col] for r in records])
    out["Diptera_relative_to_Hymenoptera"]={
        "source_sites_increasing_fly_vs_bee_visit_balance":int(np.sum(v>0)),
        "ratio_of_lit_dark_visit_ratios_geometric_equal_site":float(np.exp(v.mean())),
        "source_5or6_unit_ratio_range":list(map(float,[np.exp(v.min()),np.exp(v.max())])),
        "does_not_identify_true_flower_pollen_compensation":True}
    return out

def main(out):
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite previously completed original scientific result")
    swiss,ssource=standardize_swiss()
    tibet,tsource=standardize_tibet()
    tday=summarize(tibet["Diurnal"])
    snorm=summarize(swiss,cover=True)
    sraw=summarize(swiss)
    tnight=summarize(tibet["Nocturnal"])
    across={}
    for guild in ORDERS:
        # Difference between TWO individual experiments, not random-effects
        # meta-analysis. Uncertainty is only within-study site resampling.
        across[guild]={
          "Swiss_minus_Tibetan_equal_site_daytime_log_visitation_ratio":
              float(sraw[guild]["equal_site_log_visit_ratio"]-
                    tday[guild]["equal_site_log_visit_ratio"]),
          "Swiss_minus_Tibetan_daytime_log_ratio_after_Swiss_cover_adjustment":
              float(snorm[guild]["equal_site_log_visit_ratio"]-
                    tday[guild]["equal_site_log_visit_ratio"]),
          "two_primary_studies_only":True,
          "Swiss_2016_vs_Tibetan_2024_light_spectra_flora_geography_not_equivalent":True,
          "no_p_value_or_global_mean_misrepresented":True}
    result={"schema":"ALAN_20261114_two_original_lit_dark_field_guild_composition_flowercover_qa_v1",
       "status":"REAL_ORIGINAL_TWO_STUDIES_SITE_BLOCK_HARMONIZED_VISITS_NO_EFFECTIVE_POLLEN_META",
       "sources":{
           "Swiss_Giavi2021":ssource,"Tibet_Lu2026":tsource,
           "Swiss_original_C_Ctrl_T_Treatment_vs_Tibetan_Lit_Unlit_verified":True,
           "published_independent_original_papers":2,
           "distinct_broad_field_regions":2,
           "observed_plant_flower_cover_denominator_only_in_Swiss":True},
       "site_block_equal_results":{
         "Swiss2016_daytime_per_survey_original":sraw,
         "Swiss2016_daytime_per_recorded_plant_cover":snorm,
         "Tibet2024_daytime_per_survey_original":tday,
         "Tibet2024_nocturnal_per_survey_original":tnight},
       "cross_field_difference_NOT_meta":across,
       "source_methodological_caveats":{
           "not_global_ALAN_pollination_syndrome_causal_effect":True,
           "Swiss_plant_cover_proxy_ascertained_at_original_source_event":True,
           "Lu_no_plant_flower_availability_count":True,
           "site_block_bootstrap_not_two_study_meta_analysis":True,
           "all_guild_groups_taxonomic_orders_not_evolved_flower_pollination_syndromes":True,
           "no_direct_day_night_local_phenology_mismatch":True,
           "no_predator_attack_or_pollen_replacement_measured":True,
           "no_island_mainland_experiment":True,
           "all_source_plant_site_event_rows_excluded_from_output":True,
           "missing_visitor_guild_positives_not_absence":True}
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_TWO_ORIGINAL_FIELD_STUDY_DAY_GUILD_FLORAL_COVER",
          json.dumps({"status":result["status"],"sources":result["sources"],
                      "Swiss_site_normalized":snorm,
                      "Swiss_raw":sraw,"Tibet_day":tday,
                      "Tibet_night":tnight,"cross_site":across},
                      ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    main(a.out)
