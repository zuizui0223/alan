#!/usr/bin/env python3
"""Giavi 2021 original Zenodo 6-meadow true source plant x pair support gate.

Official Zenodo record 4540407 data.Rdata md5:f60b21d82c163cbb53520fda0bc61985.
Published source R code lmer: dVis.log1 ~ plant:treatment + plab.log.z +
(1|pair), where pair is site x sampling repeat. Keep actual numeric zeros
reported by the authors; missing plant/guild/arm not imputed to zero.
Plant × site x repeat pair co-observation with independently sourced
Kew PLANTATT and Doring moisture. Do not fit interaction if fewer than
3 original meadow sites contain TWO wet and TWO drier species among matched
identifiable named plants. A raw source row is not an independent lamp.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,tempfile,urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
import pyreadr
from scripts.analysis.alan_anic2022_13species_kew_plantatt_20261120 import (
    pull_kew_binary,trait_lookup
)
from scripts.analysis.alan_giavi2021_plantatt_moisture_dayvisits_20261120 import (
    source_name,add_doring_only_records
)
SCHEMA="alan_Giavi2021_true_Zenodo_original6meadow_pair_floral_resource_source_gate_v1"
ORIGINAL="https://zenodo.org/records/4540407/files/data.Rdata?download=1"
MD5="f60b21d82c163cbb53520fda0bc61985"
def check(ok,msg):
    if not ok:raise RuntimeError("ORIGINAL_GIAVI_ZENODO_SIXMEADOW_SOURCE_STOP: "+msg)
def fetch():
    req=urllib.request.Request(ORIGINAL,headers={
        "User-Agent":"ALAN-source-authentic-ecology-original-grognuz/2026"})
    with urllib.request.urlopen(req,timeout=45) as r:
        raw=r.read(150000)
    check(hashlib.md5(raw).hexdigest()==MD5,"official Giavi RData MD5 mismatch")
    with tempfile.NamedTemporaryFile(suffix=".Rdata") as f:
        f.write(raw);f.flush()
        obj=pyreadr.read_r(f.name)
    check("DAT" in obj,"original Giavi 'DAT' not in RData")
    return obj
def inspect(original,traits):
    f=original["DAT"].copy()
    columns=["treatment","plant","insect_group","site","repeate","plab","dVis","pair"]
    check(set(columns).issubset(f),"source field schema changed")
    for key in ("treatment","plant","insect_group","site","repeate","pair"):
        f[key]=f[key].astype(str)
    for key in ("plab","dVis"):
        f[key]=pd.to_numeric(f[key],errors="coerce")
    check(f.treatment.isin(["C","T"]).all(),"unexpected original ALAN experimental code")
    check(f["dVis"].dropna().ge(0).all(),"negative source visit counts")
    check(f["plab"].dropna().ge(0).all(),"negative source floral resources")
    g=f.groupby(["site","repeate","plant","insect_group","treatment"],dropna=False).agg(
        source_rows=("dVis","size"),source_unique_visits=("dVis","sum"),
        original_flower_mean=("plab","mean")).reset_index()
    observed=g.groupby(["site","repeate","plant","insect_group"]).treatment.nunique()
    complete=observed[observed==2].reset_index()[
        ["site","repeate","plant","insect_group"]]
    any_arm=g.groupby(["site","repeate","plant"]).treatment.nunique()
    true_arm=any_arm[any_arm==2].reset_index()[["site","repeate","plant"]]
    # Species should not be labelled as wet if missing the original UK score.
    named=set(f.plant.unique())
    source_keys={}
    for name in sorted(named):
        binom=source_name(name)
        trait=traits.get(binom) if binom else None
        source_keys[name]={
            "original_source_plant_label":name,
            "fully_resolved_binomial":binom,
            "F_PLANTATT":float(trait["F_PLANTATT"]) if trait and trait["F_PLANTATT"] is not None else None,
            "F_Doring":float(trait["F_Doring"]) if trait and trait["F_Doring"] not in (None,"NA","") else None,
            "source_original_sampled_light_arms":sorted(f.loc[f.plant.eq(name),"treatment"].unique().tolist()),
            "original_meadow_pairs_both_arms":int(true_arm[true_arm.plant.eq(name)].site.nunique()),
            "original_plant_is_not_independent_lamp":True}
    per_site={}
    for score in ("F_PLANTATT","F_Doring"):
        group=collections.defaultdict(lambda:{"wet":set(),"drier":set(),
                                              "both_arms_species":set()})
        for row in true_arm.itertuples(index=False):
            source=source_keys[row.plant]
            trait=source[score]
            if trait is None:continue
            x=group[row.site]
            x["both_arms_species"].add(row.plant)
            x["wet" if trait>7 else "drier"].add(row.plant)
        per_site[score]={
            "meadows_with_any_wet_and_drier":sum(bool(v["wet"]) and bool(v["drier"]) for v in group.values()),
            "meadows_with_2_wet_and_2_drier":sum(len(v["wet"])>=2 and len(v["drier"])>=2 for v in group.values()),
            "source_meadows_with_scored_species":len(group),
            "source_unique_scored_species_with_2plus_meadow_arms":
               sum(x[score] is not None and x["original_meadow_pairs_both_arms"]>=2
                   for x in source_keys.values())}
    return {
        "schema":SCHEMA,
        "source_original_RData_md5":MD5,
        "source_original_rows":len(f),
        "source_original_column_names":list(f),
        "source_original_distinct_meadow_pairs":int(f.site.nunique()),
        "source_original_sampling_repeats":int(f[["site","repeate"]].drop_duplicates().shape[0]),
        "source_original_plant_taxon_labels":int(f.plant.nunique()),
        "source_original_insect_groups":sorted(f.insect_group.unique().tolist()),
        "source_original_treatment_levels":sorted(f.treatment.unique().tolist()),
        "original_source_missing_visit_rows":int(f.dVis.isna().sum()),
        "original_source_zero_visit_rows":int(f.dVis.eq(0).sum()),
        "original_source_missing_flower_abundance_rows":int(f.plab.isna().sum()),
        "original_source_zero_flower_abundance_rows":int(f.plab.eq(0).sum()),
        "original_matched_site_repeat_plant_any_guild":len(true_arm),
        "original_matched_site_repeat_plant_insectguild":len(complete),
        "source_original_plant_lookup":source_keys,
        "source_meadow_moisture_support":per_site,
        "within_site_wetness_response_model_fitted":False,
        "source_original_ALAN_plant_moisture_mechanism_causal_mediation_fitted":False,
        "no_assumed_zero_for_unreported_plants":True,
        "original_flower_abundance_is_measured_and_must_be_adjusted":True,
        "source_raw_plant_visitation_records_not_committed":True,
    }
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    original=fetch()
    b=pull_kew_binary();traits,_=trait_lookup(b)
    traits=add_doring_only_records(b,traits)
    report=inspect(original,traits)
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GIAVI2021_ORIGINAL_RDATA_PAIR_FLORAL_ABUNDANCE_SOURCE_GATE "+
          json.dumps(report,ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
