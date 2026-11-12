#!/usr/bin/env python3
"""Actual species-direct floral morphology x ALAN, with source/rights guardrails.

Real input:
- frozen ALAN 10,698 original plant-photo cell years, same-year VIIRS, WCVP;
- original versus v0.11 GloBI moth-positive plant documentation labels;
- island source-scale 2026-09-08 species-DIRECT high/medium floristic traits.

No trait imputation, no inference of local moth pollination, no source rows published.
A moth-like floral form is a morphology PROXY, not effective moth pollination.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from alan_gecko_primary_fiveyear import primary
from alan_globi_crossland_light_robust_20261111 import fit

ROOT = Path(__file__).resolve().parents[2]
ID = "wcvp_accepted_plant_name_id"
BETA = "beta_island_minus_mainland_slope_pp_per_log1p"
CLASSES = ["old_only", "common", "new_only", "neither"]
FLORAL_FORMS = frozenset(["tubular", "salverform", "funnel_trumpet"])
TRAITS = ("flower_primary_color", "floral_form", "tube_depth_class")
LEDGER_SHA = "7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577"


def read_direct_traits(path: Path, high_only=False):
    if hashlib.sha256(path.read_bytes()).hexdigest() != LEDGER_SHA:
        raise ValueError("Unexpected source-scale island trait ledger digest")
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    assert (len(d), d.accepted_species.nunique()) == (134268, 62423)
    assert set(TRAITS) <= set(d.trait_name)
    d = d[d.trait_name.isin(TRAITS) & d.quality.isin(["high"] if high_only else ["high", "medium"])
          & d.resolution_status.eq("resolved")].copy()
    if d.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Ambiguous source trait per accepted name")
    # Exact single-state direct species value; missing/multistate = genuinely UNKNOWN.
    d = d[d.multistate.eq("False") & d.normalized_value.ne("")
          & ~d.normalized_value.str.contains("|", regex=False)].copy()
    pivot = d.pivot(index="accepted_species", columns="trait_name", values="normalized_value")
    for col in TRAITS:
        if col not in pivot:
            pivot[col] = pd.NA
    return pivot[list(TRAITS)]


def name_labels(sp: pd.DataFrame, old: pd.DataFrame):
    """Published old direct accepted-name match; new full-v0.11 accepted-ID match."""
    flags = pd.read_csv(ROOT / "results/globi_v011_fiveyear/source_moth_species.csv",
                        dtype={"accepted_id": str})
    assert flags.accepted_id.is_unique
    news = set(flags.loc[flags.v011_wcvp.eq(True), "accepted_id"])
    olds = set(old.loc[old.guild.eq("moth_family_visit"), "plant_name"])
    names = sp[[ID, "wcvp_accepted_name"]].drop_duplicates().copy()
    assert names[ID].is_unique
    names["old"] = names.wcvp_accepted_name.isin(olds)
    names["new"] = names[ID].astype(str).isin(news)
    names["class"] = np.select(
        [names.old & ~names.new, names.old & names.new,
         ~names.old & names.new],
        ["old_only", "common", "new_only"], default="neither")
    observed = names["class"].value_counts().to_dict()
    assert all(observed.get(k, 0) == v for k, v in
               [("old_only", 713), ("common", 612), ("new_only", 132), ("neither", 7124)])
    return names


def classify(names, trait_table):
    names = names.merge(trait_table, left_on="wcvp_accepted_name", right_index=True,
                        how="left", validate="one_to_one")
    color = names.flower_primary_color
    form = names.floral_form
    depth = names.tube_depth_class
    names["known_color"] = color.notna()
    names["known_form"] = form.notna()
    names["known_depth"] = depth.notna()
    names["eligible_broad"] = names.known_color & names.known_form
    names["eligible_strict"] = names.eligible_broad & names.known_depth
    names["white_single"] = color.eq("white").fillna(False)
    names["form_tube_single"] = form.isin(FLORAL_FORMS)
    names["broad_candidate"] = (names.eligible_broad & names.white_single &
                                names.form_tube_single)
    names["strict_candidate"] = names.broad_candidate & names.eligible_strict & depth.eq("deep").fillna(False)
    names["broad_other_known"] = names.eligible_broad & ~names.broad_candidate
    names["strict_other_known"] = names.eligible_strict & ~names.strict_candidate
    assert (names.broad_candidate & ~names.eligible_broad).sum() == 0
    assert (names.strict_candidate & ~names.eligible_strict).sum() == 0
    assert (names.eligible_broad == (names.broad_candidate | names.broad_other_known)).all()
    assert (names.eligible_strict == (names.strict_candidate | names.strict_other_known)).all()
    return names


def taxon_support(names):
    cols=["known_color", "known_form", "known_depth", "eligible_broad",
          "eligible_strict", "white_single", "form_tube_single",
          "broad_candidate", "strict_candidate", "broad_other_known"]
    ret={}
    for name in CLASSES:
        part=names.loc[names["class"].eq(name)]
        ret[name]={"species":int(len(part)),
                   **{p:int(part[p].sum()) for p in cols}}
    return ret


def stratum_test(names, target, a="old_only", b="common", gate=None):
    """Post-hoc morphology enrichment. Null = equal rates among *documented traits*."""
    subsets = []
    for tag in [a,b]:
        v = names.loc[names["class"].eq(tag)]
        if gate is not None:
            v = v.loc[v[gate]]
        subsets.append(v)
    x,y = subsets
    n1,n2=int(x[target].sum()),int(y[target].sum())
    odds,p=fisher_exact([[n1,len(x)-n1],[n2,len(y)-n2]],alternative="two-sided") if len(x) and len(y) else (np.nan,np.nan)
    return {"trait_flag":target, "scope":gate or "all_accepted_source_species",
            "old_only_positive":n1, "old_only_denominator":int(len(x)),
            "common_positive":n2, "common_denominator":int(len(y)),
            "odds_ratio":float(odds) if np.isfinite(odds) else None,
            "p_nominal_posthoc":float(p) if np.isfinite(p) else None}


def site_count(s, obs, mask):
    return (obs.loc[mask].groupby("key").size()
            .reindex(s.key, fill_value=0).to_numpy(dtype=np.int64))


def regression(s, count, name, ratio_denominator=None):
    """Same original year x geography FE, same island/mainland weights and controls.

    If ratio_denominator is not None, use only sites with >=2 independently
    classified species and both island/mainland same-year geographic support,
    adjusting logged denominator. Return descriptive post-hoc coefficient only
    if geographic/survey support is sufficient.
    """
    if ratio_denominator is None:
        denom=s.n_species.to_numpy(dtype=np.int64)
        g=s.copy()
        response=count
    else:
        if not np.all(count <= ratio_denominator):
            raise ValueError("Candidate species exceed independently classified species")
        selected=ratio_denominator >= 2
        g=s.iloc[np.flatnonzero(selected)].copy()
        response=count[selected]
        denom=ratio_denominator[selected]
        g["n_species"] = denom
        g["logspecies"] = np.log1p(denom)
    n_pos_island=int(np.sum((response>0) & g.land_class.eq("island").to_numpy()))
    n_pos_main=int(np.sum((response>0) & g.land_class.eq("mainland").to_numpy()))
    both_blocks=int(g.loc[response>0].groupby("region5").land_class.nunique().eq(2).sum())
    all_blocks=int(g.groupby("region5").land_class.nunique().eq(2).sum())
    record={
        "response":name, "denominator":"original_all_photo_species" if ratio_denominator is None else "same_site_known_morph_species_2plus",
        "positive_island_cell_years":n_pos_island,
        "positive_mainland_cell_years":n_pos_main,
        "positive_geographic_blocks_with_both_classes":both_blocks,
        "total_photo_site_years":len(g),
        "total_island_site_years":int(g.land_class.eq("island").sum()),
        "matched_geographic_blocks":all_blocks,
    }
    if not n_pos_island or not n_pos_main or len(g)<50 or all_blocks<10:
        record.update(model_status="insufficient_common_support_for_regression", beta=None)
        return record
    g["all_moth"] = response
    fit_model=fit(g,"all","moth",balance=True)
    if not fit_model.get("estimable"):
        record.update(model_status="not_estimable", fit=fit_model)
        return record
    record.update(model_status="EXPLORATORY_observed_photo_composition_only",
                  fit=fit_model,
                  biologically_interpretable=(n_pos_island>=20 and
                                              n_pos_main>=20 and both_blocks>=10))
    return record


def run(inputs:Path,signed:Path,ledger:Path,out:Path):
    if out.exists():
        raise ValueError("Refuse overwriting existing source-signed result")
    sites, sha, sp, old = primary(inputs,signed)
    assert len(sites)==10698 and int(sites.island.sum())==933 and sites.region5.nunique()==78
    assert sites.key.is_unique
    base=fit(sites, "all", "moth",balance=True)
    np.testing.assert_allclose(base[BETA],-1.9017206220556133,atol=1e-9,rtol=0)
    accepted=name_labels(sp,old)
    labels={"direct_high_medium":classify(accepted,read_direct_traits(ledger)),
            "direct_high_only":classify(accepted,read_direct_traits(ledger,high_only=True))}
    results={}
    for qual,names in labels.items():
        traits=taxon_support(names)
        obs=sp[["key",ID]].merge(names[[ID, "class", "known_color", "known_form",
           "known_depth", "eligible_broad", "eligible_strict", "white_single",
           "form_tube_single", "broad_candidate", "strict_candidate", 
           "broad_other_known", "strict_other_known"]],on=ID,validate="many_to_one")
        assert not obs.duplicated(["key",ID]).any()
        assert obs[ID].nunique()==8581
        checks=site_count(sites,obs,obs[ID].notna())
        assert np.array_equal(checks, sites.n_species.to_numpy(dtype=np.int64))
        masks={
          "white_single": obs.white_single,
          "form_tube_single":obs.form_tube_single,
          "trait_color_form_known":obs.eligible_broad,
          "trait_color_form_depth_known":obs.eligible_strict,
          "mothlike_white_tube_broad":obs.broad_candidate,
          "mothlike_white_tube_deep_strict":obs.strict_candidate,
          "other_known_color_form":obs.broad_other_known,
          "old_only_mothlike_white_tube":obs["class"].eq("old_only") & obs.broad_candidate,
          "common_mothlike_white_tube":obs["class"].eq("common") & obs.broad_candidate,
          "new_only_mothlike_white_tube":obs["class"].eq("new_only") & obs.broad_candidate,
          "neither_mothlike_white_tube":obs["class"].eq("neither") & obs.broad_candidate,
          "old_only_eligible_broad":obs["class"].eq("old_only") & obs.eligible_broad,
          "common_eligible_broad":obs["class"].eq("common") & obs.eligible_broad,
          "new_only_eligible_broad":obs["class"].eq("new_only") & obs.eligible_broad,
        }
        counts={k:site_count(sites,obs,mask) for k,mask in masks.items()}
        assert np.array_equal(counts["mothlike_white_tube_broad"],
                              sum((counts[x] for x in [
                              "old_only_mothlike_white_tube","common_mothlike_white_tube",
                              "new_only_mothlike_white_tube","neither_mothlike_white_tube"])))
        assert np.array_equal(counts["trait_color_form_known"],
                              counts["mothlike_white_tube_broad"]+
                              counts["other_known_color_form"])
        regress={k:regression(sites,v,k) for k,v in counts.items()}
        regress["broad_candidate_fraction_among_known"]=regression(
            sites,counts["mothlike_white_tube_broad"],"broad_candidate_fraction_among_known",
            counts["trait_color_form_known"])
        regress["strict_candidate_fraction_among_known"]=regression(
            sites,counts["mothlike_white_tube_deep_strict"],"strict_candidate_fraction_among_known",
            counts["trait_color_form_depth_known"])
        comparisons=[
            stratum_test(names,"eligible_broad", gate=None),
            stratum_test(names,"eligible_strict",gate=None),
            stratum_test(names,"broad_candidate",gate="eligible_broad"),
            stratum_test(names,"strict_candidate",gate="eligible_strict"),
            stratum_test(names,"white_single",gate="known_color"),
            stratum_test(names,"form_tube_single",gate="known_form"),
        ]
        results[qual]={
            "source_species_coverage_by_globi_label":traits,
            "species_old_only_vs_common_descriptive_fisher_tests":comparisons,
            "photo_site_year_numerators_and_models":regress,
            "total_photo_species_site_year_units":int(len(obs)),
        }
    result={
        "schema":"ALAN_direct_species_flower_color_form_tube_vs_source_nightlight_v1",
        "status":"REAL_SHA_VERIFIED_DIRECT_TRAIT_SOURCE_EXPLORATORY",
        "source":{
            "island_trait_snapshot":{"source_repo":"zuizui0223/island",
                "run_id":34191508045,"artifact_id":10042377593,
                "ledger_sha256":LEDGER_SHA,
                "evidence":"source-species direct high/medium, not genus/LLM proxy, not public to redistribute"},
            "original_ALAN_5year_target_sha":sha,
            "ALAN_base_old_source_beta":base[BETA],
            "WCVP_match":"exact accepted scientific species NAME only; unmatched remain unknown"},
        "study":{"photo_site_years":len(sites),
            "island_site_years":int(sites.island.sum()),
            "mainland_site_years":int((~sites.land_class.eq("island")).sum()),
            "accepted_plant_species":int(sp[ID].nunique()),
            "independent_5deg_geographic_blocks":int(sites.region5.nunique()),
            "years":[2019,2021,2023,2024,2025],
            "ALAN":"original same-year four sampled monthly VIIRS brightness, not in situ ground lux",
            "inferential_unit":"original distinct photo-accepted species per 0.05degree cell x year",
        },
        "results":results,
        "limitations":[
            "Even full white+long tube+deep trait is an unvalidated MOTH-LIKE PHENOTYPE, not effective nocturnal moth pollination",
            "Trait-not-recorded is UNKNOWN, never coded as negative",
            "Species high/medium data are independently curated and not a globally random morphological sample",
            "No proof of local moth presence/absence, plant extinction, pollinator shift or ALAN causation",
            "Source-level old-only GloBI labelling is not biological absence in newer GloBI",
            "Global source trait ledger may not be redistributed; output includes only aggregate stats and no species x site x trait rows",
            "Site-year effects are exploratory after studying historic GloBI label contrasts; p values unadjusted"
        ],
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("inputs","signed","ledger","output"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    result=run(a.inputs,a.signed,a.ledger,a.output)
    print(json.dumps({"status":result["status"],
                      "trait_coverage":{quality:data["source_species_coverage_by_globi_label"]
                      for quality,data in result["results"].items()},
                      "morphology_models":{
                        quality:{k:v for k,v in data["photo_site_year_numerators_and_models"].items()
                        if k in ["mothlike_white_tube_broad","mothlike_white_tube_deep_strict",
                         "broad_candidate_fraction_among_known","strict_candidate_fraction_among_known"]}
                      for quality,data in result["results"].items()}
                      },ensure_ascii=False,allow_nan=False))


if __name__=="__main__":
    main()
