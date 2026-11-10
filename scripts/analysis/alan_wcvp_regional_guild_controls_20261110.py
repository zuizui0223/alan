#!/usr/bin/env python3
"""On the SAME signed regional-WCVP plant cohort, compare GloBI moth, butterfly,
and nocturnal moth pollen-carriage documentation against prior PCNL light.

Input is genuine original 93,209,351-photo-ID plant×site-year source joined to
Kew TDWG3 native/introduced range categories. A category means botanical-
country species distribution, NOT island nativeness or wild photographed plant.
Moth/butterfly GloBI flags are published source-claims, not pollination function.
Walton nocturnal pollen carriage is NOT effective stigmatic pollen deposition.
All comparisons are post-hoc exploratory.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis import alan_full93m_moth_PCNLPrior_two_way_FE_20261109 as base
from scripts.analysis import alan_real_genus_prevalence_bias_20261110 as genus

SOURCE_SHA = "2f30984f627425c7dff96d7502c30bb1d28d8291835961fc059e8956204bc457"
N_SOURCE_SITES = 77647
N_SOURCE_ISLANDS = 562
SOURCE_SPECIES_SITE_RECORDS = 2012925
STATES = ("native_at_TDWG3", "introduced_at_TDWG3", "unknown")
GUILDS = ("moth", "butterfly", "night_pollen")
SCHEMA = "ALAN_regional_WCVP_moth_butterfly_nightpollen_matched_PCNLPrior_v1"


def require(ok, msg):
    if not ok: raise ValueError("ALAN_TDWG3_GUILD_NEGATIVE_CONTROL: "+str(msg))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda: f.read(4 << 20), b""):h.update(data)
    return h.hexdigest()


def load_sites(path):
    require(sha(path)==SOURCE_SHA, "Not SHA-authenticated original TDWG3 source-site table")
    frame=pd.read_csv(path)
    require(len(frame)==N_SOURCE_SITES and
            frame.island.nunique()==N_SOURCE_ISLANDS, "Original 562 island cohort drift")
    frame["island"]=frame.island.astype(str)
    require(not frame[["island","latbin","lonbin","year"]].duplicated().any(),
            "Source input duplicated original photo site-year")
    den = sum(frame["recorded_species__"+j] for j in STATES)
    moth = sum(frame["recorded_moth__"+j] for j in STATES)
    require(den.eq(frame.all_accepted_recorded_species).all() and
            moth.eq(frame.moth_visitor_positive_recorded_species).all() and
            int(den.sum())==SOURCE_SPECIES_SITE_RECORDS,
            "WCVP native/introduced/unknown source totals not conserved")
    require(frame.historic_ntl_complete.all() and
            frame.pcnl_years_observed.eq(frame.pcnl_years_expected).all() and
            frame.pcnl_years_expected.eq(frame.year-1992).all(),
            "Historical PCNL observation-year leakage")
    return frame


def source_filter(frame, status, min_species=3):
    den="recorded_species__"+status
    x=frame.loc[frame[den]>=min_species].copy()
    if x.empty:return x
    # Same factual eligibility as the regional parent PR: at least two
    # within-island geographically distinct cells in one observation year,
    # with mean prior log-PCNL difference >=0.02.
    valid=x.groupby(["island","year"],observed=True)[base.LIGHT].agg(["size","min","max"])
    keep=set(valid.loc[valid["size"].ge(2)&
            (valid["max"]-valid["min"]).ge(.02)].reset_index().island)
    return x.loc[x.island.isin(keep)].copy()


def fit_recorded_fraction(data, numerator, denominator, island_year=False):
    if len(data)<100 or data.island.nunique()<30:
        return {"status":"not_estimable","n_site_years":len(data),
                "n_islands":int(data.island.nunique())}
    x=data.copy()
    x["frac_pp"]=100*x[numerator]/x[denominator]
    x["log_obs"]=np.log1p(x.all_accepted_observation_uuids)
    x["log_species"]=np.log1p(x[denominator])
    if island_year:
        weight=1.0/x.groupby("island",observed=True).island.transform("size").to_numpy(float)
        iy=x.island.astype(str).str.cat(x.year.astype(str),sep="|").to_numpy()
        fit=genus.fixed_effects_regression(
            x[numerator].to_numpy(float)/x[denominator].to_numpy(float),
            x[base.LIGHT].to_numpy(float),x.log_obs.to_numpy(float),
            weight,x.island.to_numpy(),[iy])
        return {"status":"estimated_exploratory_observational",**fit,
                "n_site_years":len(x),"n_islands":int(x.island.nunique())}
    fit=base.fit(x,"island_equal")
    require(fit.get("status")=="estimated_exploratory_observational",
            "Parent source panel cannot fit identical photo-site FE")
    return {"status":"estimated_exploratory_observational",
            **fit["terms"][base.LIGHT],
            "n_site_years":len(x),"n_islands":int(x.island.nunique())}


def run(path,out):
    data=load_sites(path)
    results={}
    for status in STATES:
        subset=source_filter(data,status)
        den="recorded_species__"+status
        record={"n_source_sites":len(subset),
                "n_source_islands":int(subset.island.nunique()),
                "n_accepted_species_site_years":int(subset[den].sum()),
                "models":{}}
        for guild in GUILDS:
            number="recorded_"+guild+"__"+status
            record["models"][guild]=fit_recorded_fraction(subset,number,den)
        # A paired moth-minus-butterfly source-label contrast has exactly the
        # same denominator, same islands, same years and same error structure.
        subset["direct_moth_minus_butterfly"]=(
            subset["recorded_moth__"+status]-
            subset["recorded_butterfly__"+status])
        record["models"]["moth_minus_butterfly_direct"]=fit_recorded_fraction(
            subset,"direct_moth_minus_butterfly",den)
        moth=record["models"]["moth"]; butterfly=record["models"]["butterfly"]
        contrast=record["models"]["moth_minus_butterfly_direct"]
        if moth["status"]=="estimated_exploratory_observational":
            require(abs(moth["beta_pp"]-butterfly["beta_pp"]-
                       contrast["beta_pp"])<1e-8,
                    "Direct same-source GloBI guild coefficient difference failed")
        results[status]=record
    native=results["native_at_TDWG3"]["models"]["moth"]
    introduced=results["introduced_at_TDWG3"]["models"]["moth"]
    require(abs(native["beta_pp"]-0.7721703383470593)<1e-6 and
            abs(introduced["beta_pp"]-(-0.32651954703074604))<1e-6,
            "Parent PR28 WCVP regional source estimates not numerically reproduced")
    # Hold island×year source exposure and all-photo WCVP source denominator
    # fixed, to test which regional-status recorded flora share covaries with
    # nightlight rather than asserting cultivated plant status.
    composition={}
    for status in STATES:
        frame=data.copy()
        numerator="recorded_species__"+status
        composition[status]=fit_recorded_fraction(
            frame,numerator,"all_accepted_recorded_species")
    require(abs(sum(x["beta_pp"] for x in composition.values()))<1e-8,
            "Regional WCVP share source fractions do not sum to one")
    # Same-site native-minus-introduced test: both groups have >=3 documented
    # accepted species; do not compare disjoint regional model p values.
    both=data.loc[data.recorded_species__native_at_TDWG3.ge(3)&
                  data.recorded_species__introduced_at_TDWG3.ge(3)].copy()
    valid=both.groupby(["island","year"],observed=True)[base.LIGHT].agg(["size","min","max"])
    keep=set(valid.loc[valid["size"].ge(2)&
             (valid["max"]-valid["min"]).ge(.02)].reset_index().island)
    both=both.loc[both.island.isin(keep)].copy()
    paired={"n_source_sites":len(both),"n_source_islands":int(both.island.nunique()),
            "models":{}}
    for guild in GUILDS:
        # direct fraction difference, not divided by a second denominator
        both["response_fraction_native_minus_introduced"]=(
            both["recorded_"+guild+"__native_at_TDWG3"]/
            both.recorded_species__native_at_TDWG3-
            both["recorded_"+guild+"__introduced_at_TDWG3"]/
            both.recorded_species__introduced_at_TDWG3)
        both["one"]=1
        paired["models"][guild]=fit_recorded_fraction(
            both,"response_fraction_native_minus_introduced","one")
    result={"schema":SCHEMA,"source_SHA256":SOURCE_SHA,
            "source_original_photo_IDs":93209351,
            "source_site_years":len(data),
            "source_islands":int(data.island.nunique()),
            "regional_guilds":results,
            "native_minus_introduced_identical_sources":paired,
            "regional_composition_share_vs_PCNL":composition,
            "GloBI_positive_claim_not_true_moth_dependence":True,
            "night_pollen_carriage_not_effective_pollen_deposition":True,
            "WCVP_regional_native_not_local_island_native":True,
            "cultivated_or_wild_individual_unknown":True,
            "all_additional_models_exploratory_posthoc":True}
    folder=Path(out)
    require(not folder.exists(),"Do not replace existing source audit")
    folder.mkdir(parents=True)
    (folder/"ALAN_regional_native_GloBI_guild_controls_20261110.json").write_text(
        json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_REGIONAL_WCVP_MOTH_BUTTERFLY_NOCTURNAL_SOURCE_CONTROLS",
          json.dumps(result,sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.sites,args.out)


if __name__=="__main__":main()
