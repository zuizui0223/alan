#!/usr/bin/env python3
"""Giavi 2021 Swiss ALAN independent original plant × daytime guild source model.

Original open data Zenodo 4540407, Rdata original MD5 verified.
This completes a second independent geographical field manipulation study
for plant x insect-order composition with different site-year protocols than
Lu 2026 Tibetan study.

Unit: six original paired sites and 37 paired survey events, NOT 1884
independent treatments. Counts are in published Rdata object 'DAT';
plant opportunity is observed 'plab', not actual floral nectar supply.
Efficacy/seed and timing-mismatch are NOT measured.

Only aggregate by site (not source site IDs), taxonomic order and
source-positive plant opportunity counts. Never release the original
plant x site x visit-event rows into this repository.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import collections
import io
import itertools
import json
import math
import tempfile

import numpy as np
import pandas as pd
import pyreadr

from scripts.analysis.alan_Giavi2021_Swiss_original_data_source_audit_20261114 import (
    load_file, SRC, dataframe_info,
)

KEY=["site","pair","plant","insect_group","treatment"]
GUILDS=("Coleoptera","Diptera","Hymenoptera")
SOURCE_TOTAL=1884


def require(test, msg):
    if not bool(test):
        raise ValueError("ALAN_Switzerland_Giavi2021_ORIGINAL_SOURCE: "+msg)


def full_six_block_signflip(values):
    v=np.asarray(values,float)
    require(len(v)==6 and np.isfinite(v).all(),
            "Six true independent original experimental site pair contrasts required")
    obs=abs(float(v.mean()))
    signs=itertools.product([-1.,1.],repeat=6)
    null=np.array([abs(np.mean(v*np.array(s))) for s in signs])
    return {"site_pairs":6,
      "mean_paired_source_site_contrast":float(v.mean()),
      "positive_site_differences":int(np.sum(v>0)),
      "negative_site_differences":int(np.sum(v<0)),
      "two_sided_full_64_site_signflip_fraction_descriptive":float(
           np.mean(null>=obs-1e-12)),
      "not_claimed_true_assignment_randomization_test":True}


def load_original_data():
    b,receipt=load_file("data.Rdata")
    require(b is not None and receipt["source_authenticated"],
            "Original 2016 published Rdata not retrievable or hash mismatch")
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/"data.Rdata"
        p.write_bytes(b)
        obj=pyreadr.read_r(str(p))
    require(set(obj)>={"DAT","gradient"},"Original source R objects missing")
    x=obj["DAT"].copy()
    require(len(x)==SOURCE_TOTAL and
            {"treatment","plant","insect_group","site","repeate",
             "plab","dVis","pair"}.issubset(x.columns),
            "Original R source row count or schema changed")
    return x,receipt


def check_and_prepare(a):
    a=a.copy()
    for k in ("treatment","plant","insect_group","site","pair"):
        a[k]=a[k].astype(str)
    a["dVis"]=pd.to_numeric(a.dVis,errors="raise")
    a["plab"]=pd.to_numeric(a.plab,errors="raise")
    a["repeate"]=pd.to_numeric(a.repeate,errors="raise")
    require(a[["dVis","plab","repeate"]].notna().all().all(),
            "Original measured quantitative source missing")
    require(a.dVis.ge(0).all() and a.plab.ge(0).all(),
            "Negative measured source visitation or plant abundance")
    require(set(a.treatment)=={"C","T"} and
            set(a.insect_group)==set(GUILDS),
            "Original dark C, lit T or insect guild source drift")
    require((a.plant.nunique(),a.site.nunique(),a.pair.nunique())==(21,6,37),
            "Original 21 taxa 6 sites 37 paired survey events changed")
    require(not a.duplicated(KEY).any(),
            "Source paper×plant×order×light key duplicated; do not pseudo-match")
    by_pair=a.groupby("pair").site.nunique()
    require(by_pair.eq(1).all(),"Original pair spans multiple published sites")
    both=a.groupby("pair").treatment.nunique()
    require(both.eq(2).all(),
            "Source paired event missing lit or dark condition")
    return a


def by_guild(a):
    """Within each true paired site, compare each treatment guild, no imputed visitor absence."""
    totals=a.groupby(["insect_group","treatment"],observed=True).dVis.sum().unstack("treatment")
    require(totals.index.nunique()==3 and set(totals.columns)=={"C","T"},
            "Three guilds or treatments altered")
    pairs=a.groupby(["site","pair","insect_group","treatment"],observed=True).dVis.sum().unstack("treatment")
    # All 37 paired original events × 3 recorded orders need to be represented
    # by original controls AND treatment. Missing is a DATA GAP not a zero.
    paired=pairs.loc[pairs.notna().all(axis=1)].copy()
    bysite=paired.groupby(["site","insect_group"],observed=True)[["C","T"]].sum()
    # Sites may have different numbers of the 37 paired events.
    source_design={
        "original_source_rows":len(a),
        "total_original_visitor_count":float(a.dVis.sum()),
        "original_original_field_site_pairs":int(a.site.nunique()),
        "original_paired_transect_sampling_event_count":int(a.pair.nunique()),
        "records_with_nonzero_visits":int(a.dVis.gt(0).sum()),
        "source_treatment_rows":{k:int(v) for k,v in a.treatment.value_counts().items()},
        "original_sample_pair_site_counts":{
            "min":int(a.groupby("site").pair.nunique().min()),
            "max":int(a.groupby("site").pair.nunique().max())},
        "total_paired_insect_group_event_rows":int(len(paired)),
        "missing_any_treatment_pair_guild_event_not_imputed":int(pairs.isna().any(axis=1).sum()),
    }
    groups={}
    for guild in GUILDS:
        b=bysite.xs(guild,level="insect_group")
        require(len(b)==6,"Original independent Swiss site pair missing for "+guild)
        dark,lit=[float(totals.loc[guild,k]) for k in ("C","T")]
        site_ratio=(b.T / b.C.replace(0,np.nan)) if False else None
        site_d=(b["T"]-b["C"]).to_numpy(float)
        # Site pair diagnostics even if each site sampled 6 or 7 surveys.
        # Because C and T within a pair have equal observation opportunity,
        # an aggregated site ratio is only an OBSERVED visitation count ratio.
        stat=full_six_block_signflip(site_d)
        groups[guild]={
          "original_dark_visits":dark,
          "original_lit_visits":lit,
          "site_pair_original_unlit_visits":[float(x) for x in b["C"]],
          "site_pair_original_lit_visits":[float(x) for x in b["T"]],
          "original_lit_minus_dark_visits":lit-dark,
          "observed_percent_light_change":float((lit/dark-1)*100) if dark>0 else None,
          "source_visit_difference_site_signs":stat,
          "no_effective_pollen_seed_outcome_measured":True}
    # Cannot claim the original 2384 visits if original Rdata is an analysis
    # subset that does not perfectly recreate the reported raw total.
    return source_design,groups


def original_plant_abundance_bias(a):
    # plant abundance repeated across insect_group rows; deduplicate at
    # original plant×paired sampling event×treatment before computing opportunity.
    q=a.groupby(["pair","site","plant","treatment"],observed=True)
    ranges=q.plab.agg(["min","max"])
    inconsistent=int((ranges["max"]-ranges["min"]>1e-8).sum())
    plants=q.plab.median().reset_index()
    by=plants.groupby("treatment").plab.agg(["sum","mean","median"])
    return {"source_original_plant_abundance_values":len(plants),
            "multi_guild_same_plant_event_inconsistent_abundance_values":inconsistent,
            "source_cumulative_deduplicated_plant_abundance":{
                t:float(by.loc[t,"sum"]) for t in ("C","T")},
            "number_plants_positive_in_both_treatments":int(
                (plants.groupby(["plant","treatment"]).plab.sum().unstack("treatment").fillna(0)>0).all(axis=1).sum()),
            "plant_abundance_used_as_original_author_covariate_not_true_flower_count":True}


def eligible_plant_contrast(a):
    """Only plants observed in both light conditions in at least two source site pairs."""
    visits=a.groupby(["site","pair","plant","insect_group","treatment"],observed=True).dVis.sum().unstack("treatment")
    both=visits.loc[visits.notna().all(axis=1)]
    source_support=both.reset_index().groupby("plant").site.nunique()
    supported=source_support.loc[source_support.ge(2)]
    return {"unique_original_plants":int(a.plant.nunique()),
       "plants_with_pairwise_lit_dark_positive_record_both_and_2plus_sites":
           int(len(supported)),
       "plants_with_5plus_sites_positive_source_support":
           int((source_support>=5).sum()),
       "original_plant_by_site_by_guild_within_pair_matched_counts":
           int(len(both)),
       "source_rows_undocumented_not_falsely_zero":True,
       "no_taxon_specific_effective_pollen_or_prior_moth_shift_inferred":True}


def run(out):
    out=Path(out)
    require(not out.exists(),"Original audit output already exists")
    raw,receipt=load_original_data()
    x=check_and_prepare(raw)
    source,guilds=by_guild(x)
    plant=eligible_plant_contrast(x)
    plab=original_plant_abundance_bias(x)
    result={"schema":"ALAN_Swiss_Giavi2021_6site_37pair_raw_guild_light_source_comparison_v1",
       "status":"ORIGINAL_2016_SWISS_SIX_PAIRED_FIELD_SITE_SOURCE_VERIFIED",
       "citation":{"article_DOI":"10.1038/s41467-021-22011-8",
                   "original_Zenodo_DOI":"10.5281/zenodo.4540407",
                   "original_Rdata_SHA256":receipt["sha256"],
                   "original_Rdata_MD5":receipt["original_md5"]},
       "source_design":source,
       "original_guild_visit_contrasts":guilds,
       "original_plant_opportunity":plab,
       "original_source_plant_support":plant,
       "interpretation":{
           "independent_Swiss_from_Tibetan_Lu2026_original":True,
           "same_21_taxa_not_cross_historical_shift_proof":True,
           "no_same_site_real_nocturnal_visits":True,
           "no_same_plant_fruit_seed_reproduction":True,
           "no_light_phenology_or_island_mainland_contrast":True,
           "no_actual_same_flower_predator_attack":True,
           "guild_light_effect_is_daytime_only":True,
           "no_unrecorded_flower_or_visitor_inferred_absent":True,
           "no_raw_source_plant_site_pair_rows_exported":True,
           "exact_signflip_is_DESCRIPTIVE_not_design_randomization":True,
       }}
    out.parent.mkdir(exist_ok=True,parents=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_GIAVI2021_REAL_SWISS_DAYTIME_GUILD_REANALYSIS",
          json.dumps(result,ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    run(args.out)
