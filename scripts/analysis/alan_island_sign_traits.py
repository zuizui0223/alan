"""Describe directly scored traits on both sides of the frozen island slope sign."""
import argparse
import hashlib
from pathlib import Path

import pandas as pd

from alan_globi_v011_same_fiveyear import taxonomy

ROOT = Path(__file__).resolve().parents[2] / "results"
OUT = ROOT / "species_interpretation"


def run(inputs, direct, mating):
    s = pd.read_csv(OUT / "shared_544.csv").set_index("plant")
    updates = pd.read_csv(ROOT / "mechanism_exemplars/transport_updates.csv")
    s.loc[updates.plant, "transport_review_status"] = "documented"
    assert hashlib.sha256(direct.read_bytes()).hexdigest() == "7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577"
    d = pd.read_csv(direct, dtype=str, keep_default_na=False)
    d = d[d.accepted_species.isin(s.index) & d.quality.isin(["medium","high"]) &
          d.resolution_status.eq("resolved") & d.multistate.eq("False") &
          d.normalized_value.ne("") & ~d.normalized_value.str.contains("|",regex=False)]
    axes = d.pivot(index="accepted_species", columns="trait_name",values="normalized_value").reindex(s.index)
    broad = axes.flower_primary_color.eq("white") & axes.floral_form.isin(["tubular","salverform","funnel_trumpet"])
    known = axes.flower_primary_color.notna() & axes.floral_form.notna()
    axes["white_tubular"] = broad.where(known)
    axes["white_tubular_deep"] = (broad & axes.tube_depth_class.eq("deep")).where(known & axes.tube_depth_class.notna())
    assert axes.white_tubular.count()==103 and axes.white_tubular.sum()==4 and axes.white_tubular_deep.count()==34
    assert hashlib.md5(mating.read_bytes()).hexdigest()=="c6b816501a98d336b6fa43fa60ae7ca4"
    m = pd.read_csv(mating, encoding="latin1")
    names = m.genus.str.strip() + " " + m.species.str.strip()
    cross = taxonomy(inputs,set(names.dropna()),set())
    m["plant"] = names.where(names.isin(s.index),names.map(cross.loc[cross.resolution.eq("resolved"),"accepted_name"]))
    tm = m["mean.tm"]
    m["selfing"] = (1-tm).where(tm.between(0,1))
    s["selfing_median"] = m.groupby("plant").selfing.median()
    s["genetic_tm_records"] = m.groupby("plant")["mean.tm"].count().reindex(s.index,fill_value=0)
    rows = []
    for sign,mask in [("negative",s.island_pp.lt(0)),("positive",s.island_pp.gt(0))]:
        c = s[mask]
        features = axes[mask].copy()
        features["source_compatibility"] = c.sc_source_class.replace("unknown",pd.NA)
        features["transport_documented"] = c.transport_review_status.eq("documented")
        features["PL_positive"] = c.glopl_median_log_response_ratio.gt(0).where(c.glopl_cases.gt(0))
        features["night_exclusion_loss_positive"] = c.day_only_relative_loss.gt(0).where(c.day_only_relative_loss.notna())
        for name,series in features.items():
            for value,n in series.value_counts().items():
                rows.append([sign,len(c),name,value,series.count(),n])
        for name in ["glopl_median_log_response_ratio","selfing_median"]:
            series = c[name].dropna()
            print(sign,name,len(series),series.median())
        print(sign,"any genetic tm",c.genetic_tm_records.gt(0).sum())
    result = pd.DataFrame(rows,columns=["island_sign","cohort","trait","state","scored","species"])
    result.to_csv(OUT / "island_sign_trait_counts.csv", index=False)
    cols = ["island_pp","mainland_pp","gap_pp","island_site_years","mainland_site_years","v011",
            "transport_review_status","transport_evidence","transport_source","day_only_relative_loss","night_exclusion_studies",
            "glopl_cases","glopl_median_log_response_ratio","sc_source_class","genetic_tm_records","selfing_median"]
    s.loc[updates.plant,"transport_evidence"] = updates.set_index("plant").evidence
    s.loc[updates.plant,"transport_source"] = updates.set_index("plant").study_doi
    s[cols].sort_values("island_pp").to_csv(OUT / "island_sign_reproductive_evidence.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["inputs","direct","mating"]:
        parser.add_argument(arg,type=Path)
    args = parser.parse_args()
    run(args.inputs,args.direct,args.mating)
