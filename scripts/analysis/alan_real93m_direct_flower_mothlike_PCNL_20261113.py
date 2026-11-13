#!/usr/bin/env python3
"""Independent species-direct moth-LIKE morphology × historic PCNL on 562 genuine islands.

Morphology proxy: white single-state corolla + tubular/salverform/funnel_trumpet.
Strict additionally requires documented deep tube. Not proof of moth pollination.
Retains original 93,209,351-photo ID-derived site and species populations and
source-year-correct Li-Cao 1992–2024 satellite light. Unknown traits stay UNKNOWN.

This is a post-hoc observational test, not plant change through time or causality.
Never publish per-species private trait values or per-site/private source rows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, LIGHT, load_sites, load_labels, original_accepted_species,
    site_outcome_model, species_model, require
)
from scripts.analysis.alan_direct_flower_syndrome_light_20261112 import (
    read_direct_traits, FLORAL_FORMS
)

SOURCE_SITE_YEARS = 77647
SOURCE_ISLANDS = 562
SOURCE_SPECIES_SITE_YEARS = 2012925
SOURCE_SPECIES = 14330
GROUPS = ("old_only","both","new_only","neither")


def species_morphology(species: pd.DataFrame, trait_table: pd.DataFrame,
                       v011_flags: pd.DataFrame):
    """Accepted-name crosswalk only. Missing/multistate source stays UNKNOWN."""
    taxa=species[["species", "genus", "moth_original"]].drop_duplicates("species").copy()
    require(taxa.species.is_unique and len(taxa)==SOURCE_SPECIES, "Wrong original accepted taxon set")
    for trait in ("flower_primary_color","floral_form","tube_depth_class"):
        taxa[trait]=taxa.species.map(trait_table[trait])
    taxa["known_color"]=taxa.flower_primary_color.notna()
    taxa["known_form"]=taxa.floral_form.notna()
    taxa["known_depth"]=taxa.tube_depth_class.notna()
    taxa["known_pair"]=taxa.known_color & taxa.known_form
    taxa["known_triplet"]=taxa.known_pair & taxa.known_depth
    taxa["white"]=taxa.flower_primary_color.eq("white").fillna(False)
    taxa["tubular_like"]=taxa.floral_form.isin(FLORAL_FORMS)
    taxa["broad"]=taxa.known_pair & taxa.white & taxa.tubular_like
    taxa["strict"]=taxa.known_triplet & taxa.broad & taxa.tube_depth_class.eq("deep").fillna(False)
    taxa["broad_other_known"]=taxa.known_pair & ~taxa.broad
    require(taxa.strict.le(taxa.broad).all(), "Strict must be nested inside broad")
    require((taxa.broad & ~taxa.known_pair).sum()==0,"Missing morphology wrongly positive")
    require((taxa.known_pair == (taxa.broad | taxa.broad_other_known)).all(),
            "Known classified morphology fractions not closed")
    ids = v011_flags[["accepted_name","v011_wcvp"]].drop_duplicates("accepted_name")
    require(ids.accepted_name.is_unique, "Duplicated v0.11 species accepted names")
    taxa=taxa.merge(ids,left_on="species",right_on="accepted_name",how="left",
                    validate="one_to_one")
    taxa["new"]=taxa.v011_wcvp.fillna(False).astype(bool)
    taxa["old"]=taxa.moth_original.astype(bool)
    taxa["source_class"]=np.select(
        [taxa.old & ~taxa.new, taxa.old & taxa.new, ~taxa.old & taxa.new],
        ["old_only","both","new_only"],default="neither")
    require(sum(taxa["source_class"].eq(x).sum() for x in GROUPS)==len(taxa),
            "GloBI documentation classes nonconserved")
    return taxa


def subset_support(d: pd.DataFrame, sites: pd.DataFrame, response: str):
    """Distinct islands and survey support for recorded positive plant morphology."""
    positive=d.loc[d[response].astype(bool), KEY].drop_duplicates()
    where=sites.merge(positive,on=KEY,how="inner",validate="one_to_one")
    iy=where.groupby(["island","year"],observed=True).size()
    return {
        "species_site_years_positive":int(d[response].sum()),
        "distinct_positive_photo_site_years":len(where),
        "independent_positive_islands":int(where.island.nunique()),
        "positive_island_years":int(len(iy)),
        "positive_island_years_2plus_photo_sites":int(iy.ge(2).sum()),
    }


def source_site_model(d: pd.DataFrame,sites:pd.DataFrame,attribute:str,
                      gate_islands=25, gate_sites=150):
    """Model full source cell-year denominator; no unknowns called negatives."""
    supp=subset_support(d,sites,attribute)
    if supp["independent_positive_islands"]<gate_islands or \
       supp["distinct_positive_photo_site_years"]<gate_sites:
        return {"status":"not_estimated_insufficient_independent_positive_islands",
                "positive_support":supp,
                "no_missing_traits_as_known_negative":True}
    numerator=d.groupby(KEY,observed=True)[attribute].sum().reset_index()
    matched=sites.merge(numerator,on=KEY,how="left",validate="one_to_one")
    require(matched[attribute].notna().all(),"Source original photo site lost")
    require(matched[attribute].le(matched.all_accepted_recorded_species).all(),
            "Source numerator exceeds recorded species denominator")
    out=site_outcome_model(matched,matched[attribute].to_numpy(float))
    out["status"]="exploratory_photo_documentation_composition"
    out["positive_support"]=supp
    out["response_description"]=attribute+" among ALL photo-observed accepted species, unknown true morphology retained in denominator"
    return out


def known_site_comparability(d:pd.DataFrame,sites:pd.DataFrame, key_col="known_pair"):
    """Use only site-years with >=2 distinct source-documented shape+color plants."""
    counts=d.loc[d[key_col]].groupby(KEY,observed=True).size().reset_index(name="known_n")
    tmp=sites.merge(counts,on=KEY,how="left",validate="one_to_one")
    selected=tmp.loc[tmp.known_n.ge(2)].copy()
    if selected.empty:return selected
    g=selected.groupby(["island","year"],observed=True)[LIGHT].agg(["size","min","max"])
    keep=g.loc[g["size"].ge(2)&(g["max"]-g["min"]).ge(.02)]
    selected=selected.merge(keep[[]].reset_index(),on=["island","year"],how="inner",
                            validate="many_to_one")
    return selected


def species_direct_known_only_model(d:pd.DataFrame,sites:pd.DataFrame,
                                    *, genus_fe:bool, high_only:bool=False,
                                    mixed_only:bool=False):
    """Fit only known pairwise traits: morphological noncandidate must be known.

    Genus FE tests within-genus differences, a conservative control for
    generalized urban/horticultural and taxonomic photographic preference.
    """
    selected=known_site_comparability(d,sites)
    source=d.loc[d.known_pair].copy()
    if mixed_only:
        g=source.drop_duplicates("species").groupby("genus",observed=True).broad.agg(
            ["size","sum"])
        keep=g.loc[g["size"].ge(2) & g["sum"].ge(1) &
                   (g["size"]-g["sum"]).ge(1)].index
        source=source[source.genus.isin(keep)].copy()
        selected=known_site_comparability(source,sites)
    sub=source.merge(selected[KEY],on=KEY,how="inner",validate="many_to_one")
    support={
      "known_color_form_species":int(source.species.nunique()),
      "documented_broad_mothlike_species":int(source.loc[source.broad,"species"].nunique()),
      "survey_photo_site_years_after_min2_known_and_light":int(len(selected)),
      "independent_islands_with_eligible_photo_sites":int(selected.island.nunique()),
      "independent_islands_with_mothlike_photo_sites":int(sub.loc[sub.broad,"island"].nunique()),
      "positive_photo_site_years":int(sub.loc[sub.broad,KEY].drop_duplicates().shape[0]),
      "observed_mixed_genera":int(sub.genus.nunique()),
    }
    if support["independent_islands_with_eligible_photo_sites"]<30 or \
       support["independent_islands_with_mothlike_photo_sites"]<20 or \
       support["positive_photo_site_years"]<100 or len(sub)<1000:
        return {"status":"not_estimated_insufficient_common_morphology_support",
                "source_support":support,"unknown_trait_not_implicitly_negative":True}
    try:
        fit=species_model(sub,sites,response="broad",genus_fe=genus_fe,
                          weight="island_equal")
    except ValueError as e:
        return {"status":"not_estimated_design_or_solver_failure","error":str(e),
                "source_support":support}
    fit.update(status="exploratory_direct_morphology_association",
               source_support=support,
               denominator="ONLY observed known color and form species, candidate vs measured noncandidate",
               excluded_unmeasured_species_as_missing=True,
               morphology_not_observed_effective_nocturnal_pollination=True)
    return fit


def run(sites_file:Path,globisource:Path,groups:Path,ledger:Path,
        v011:Path,out:Path)->dict:
    if out.exists():raise FileExistsError(out)
    sites=load_sites(sites_file)
    require(len(sites)==SOURCE_SITE_YEARS and sites.island.nunique()==SOURCE_ISLANDS,
            "Island and PCNL fixed source cohort changed")
    old=load_labels(globisource)
    d=original_accepted_species(groups,sites,old)
    baseline=site_outcome_model(sites, d.groupby(KEY,observed=True).moth_original.sum()
                           .reindex(pd.MultiIndex.from_frame(sites[KEY])).to_numpy(float))
    original_beta=float(baseline["terms"][LIGHT]["beta_pp"]) if "terms" in baseline else float(baseline["beta_pp"])
    require(abs(original_beta-.45722827881259)<1.e-5,
            "Fixed source original GloBI coefficient not reconstructed")

    vf=pd.read_csv(v011,dtype={"accepted_id":str})
    require(vf.accepted_name.is_unique,"Full GloBI v0.11 WCVP name source changed")

    results={}
    species_classes={}
    for quality in ("high_medium","high_only"):
        trait=read_direct_traits(ledger,high_only=quality=="high_only")
        taxa=species_morphology(d,trait,vf)
        for cls in GROUPS:
            s=taxa.loc[taxa.source_class.eq(cls)]
            species_classes[cls]=int((taxa.source_class==cls).sum()) if quality=="high_medium" else species_classes[cls]
        source=d.merge(taxa[["species","known_color","known_form","known_depth",
                    "known_pair","known_triplet","white","tubular_like",
                    "broad","strict","broad_other_known","source_class"]],
                       on="species",how="left",validate="many_to_one")
        require(len(source)==SOURCE_SPECIES_SITE_YEARS, "Species/source join lost originals")
        require(source.known_pair.notna().all(), "Source morphology unknown not classified")

        byclass={}
        for cls in GROUPS:
            a=taxa.loc[taxa.source_class.eq(cls)]
            byclass[cls]=dict(
                species=int(len(a)),
                species_known_color=int(a.known_color.sum()),
                species_known_form=int(a.known_form.sum()),
                species_known_color_form=int(a.known_pair.sum()),
                species_known_color_form_depth=int(a.known_triplet.sum()),
                species_broad_mothlike=int(a.broad.sum()),
                species_strict_mothlike=int(a.strict.sum()),
            )

        models={}
        for attr in ("known_pair","known_color","known_form","white",
                     "tubular_like","broad","strict"):
            models[attr]=source_site_model(source,sites,attr)
        # Numerator and denominator are from independent source trait records.
        models["known_only_genus_FE"]=species_direct_known_only_model(source,sites,genus_fe=True)
        models["known_only_no_genus_FE"]=species_direct_known_only_model(source,sites,genus_fe=False)
        models["known_only_mixed_genera"]=species_direct_known_only_model(
            source,sites,genus_fe=True,mixed_only=True)
        for cls in ("old_only","both","new_only"):
            attr="broad_"+cls
            source[attr]=source.broad & source.source_class.eq(cls)
            models[attr]=source_site_model(source,sites,attr,gate_islands=25,gate_sites=150)

        results[quality]=dict(
            traits_in_external_ledger=int(len(trait)),
            source_species_with_direct_color_form=int(taxa.known_pair.sum()),
            source_species_with_direct_broad_mothlike=int(taxa.broad.sum()),
            source_species_with_direct_strict_mothlike=int(taxa.strict.sum()),
            source_species_by_old_new_moth_GloBI_class=byclass,
            outcome_models=models,
            source_photo_species_years=int(len(source)),
        )
    result={
        "schema":"ALAN_93m_photo_species_direct_morphology_prior33y_PCNL_global562_v1",
        "status":"REAL_SOURCE_REPLAY_OR_STRICT_SUPPORT_FAIL_CLOSED",
        "source":{"original_photo_fullindex_n":93209351,
                  "canonical_supported_source_sites":85318,
                  "matched_islands":SOURCE_ISLANDS,
                  "matched_site_years":SOURCE_SITE_YEARS,
                  "original_WCVP_accepted_species":SOURCE_SPECIES,
                  "species_site_years":SOURCE_SPECIES_SITE_YEARS,
                  "historical_light":"Li-Cao corrected PCNL 1992–2024, arithmetic mean of yearly log(1+prior_PCNL), previous years only, original MD5 sources",
                  "same_frame_original_GloBI_beta":original_beta,
                  "direct_trait_ledger":"island run 34191508045 source-DIRECT high/medium checksum-gated, no redistributing raw trait",
                  "GloBI_source":"old cache is an observational *documentation* flag; full 0.11 change tracked separately"},
        "results":results,
        "interpretation_limits":{
           "no_photo_pixels_in_this_study":True,
           "GloBI_moth_visit_is_not_local_effective_pollination":True,
           "morphology_is_not_actual_moth_pollination":True,
           "unknown_traits_remain_unknown":True,
           "recorded_plant_proportions_not_population_change":True,
           "source_geographic_cultivated_native_bias_unresolved":True,
           "posthoc_overlap_and_multiple_comparisons":True,
           "different_ALAN_products_not_interchangeable":True
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False,default=str)+"\n")
    print("REAL_562_ISLAND_DIRECT_TRAIT_SUMMARY",
          json.dumps({q:{"known_pair":a["source_species_with_direct_color_form"],
                         "broad":a["source_species_with_direct_broad_mothlike"],
                         "strict":a["source_species_with_direct_strict_mothlike"],
                         "models":{k:{j:v.get(j) for j in ("status","beta_pp","95CI_pp","nominal_p")}
                              for k,v in a["outcome_models"].items()}}
                      for q,a in results.items()},ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("sites","globi","groups","ledger","v011","out"):
        ap.add_argument("--"+name,required=True,type=Path)
    args=ap.parse_args()
    run(args.sites,args.globi,args.groups,args.ledger,args.v011,args.out)
