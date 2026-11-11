"""Relate experimentally measured night-exclusion losses to frozen land contrasts."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from alan_gecko_primary_fiveyear import primary
from alan_globi_crossland_light_robust_20261111 import fit
from alan_globi_v011_same_fiveyear import taxonomy

ROOT = Path(__file__).resolve().parents[2]
GAP = "beta_island_minus_mainland_slope_pp_per_log1p"
METRIC = "treatment_effectiveness_metric"
ARMS = {"open pollination": "open_value", "day pollination": "day_value"}


def paired_experiments(folder, inputs, cohort):
    raw = pd.read_csv(folder / "data/DvN_LiteratureDataClean_2024-09-13.csv")
    names = set(raw.plant_species.dropna())
    tax = taxonomy(inputs, names, set())
    mapping = tax.loc[tax.resolution.eq("resolved"), "accepted_name"].to_dict()
    mapping.update({n: n for n in names & cohort})
    # Two North American studies explicitly describe white campion under this synonym.
    mapping["Silene alba"] = "Silene latifolia"
    raw["plant"] = raw.plant_species.map(mapping)
    screened = raw[raw.plant.isin(cohort) & raw.exclude.eq("n")]
    # Prespecified endpoint per study: fruit set per flower when recoverable,
    # otherwise seed production. Counts of fruits per stalk alone are omitted.
    choice = {
        "Morse_1983": "fruit set", "Anderson_2015": "seed set",
        "Bartelmess_2006": "fruit set", "Guitan_1993": "fruit set",
        "Macgregor_2019": "fruit set", "Scopece_2018": "seed set",
        "Vanderplank_2020": "fruit set", "Young_2002": "seed set",
        "Robertson_2021": "fruit set", "Cutler_2012": "fruit set",
        "Jabbari_2013": "fruit set", "Stone_2018": "fruit set",
        "Noone_2022": "fruit set", "Wonderlin_2023": "seed set",
    }
    r = screened[screened[METRIC].isin(["fruit set", "seed set"])
                 & screened.treatment_condition.isin(ARMS)].copy()
    r["arm"] = r.treatment_condition.map(ARMS)
    r["source"] = r.DOI.replace({"NA - have PDF": "10.1111/een.13296"})
    r.loc[r.study_ID.eq("Guitan_1993"), "source"] = "Guitian, Guitian & Navarro 1993, Acta Oecologica 14:219-227, Fig. 2"
    keys = ["plant", "plant_species", "study_ID", "locality", "Site_period", "year_start",
            "plant_cultivar", "additional_treatment", METRIC]
    assert not r.duplicated(keys+["arm"]).any()
    day = r[r.arm.eq("day_value")].rename(columns={"effectiveness_value": "day_value"})
    op = r[r.arm.eq("open_value")].rename(columns={"effectiveness_value": "open_value"})
    pairs = day.merge(op[keys+["open_value"]], on=keys, validate="one_to_one")
    pairs = pairs[pairs.study_ID.ne("Allison_2022")].copy()  # camera observations
    pairs["used"] = pairs[METRIC].eq(pairs.study_ID.map(choice))
    pairs["basis"] = pairs.extraction_type
    pairs.loc[pairs.study_ID.eq("Bertin_1980"), "basis"] = "Published extraction; unit unresolved, not used"
    pairs.loc[pairs.study_ID.eq("Jennersten_1991"), "basis"] = "Pods per stalk, not per flower; not used"
    assert set(pairs.loc[pairs.used, "study_ID"]) == set(choice)
    return pairs[keys+["source","basis","day_value","open_value","used"]], screened


def recent(folder):
    fruit = pd.read_csv(folder/"raspfruit.csv", sep=";").query("Plant_sp == 'RB'").groupby("Treatment").sum(numeric_only=True)
    fruit_rate = fruit.Fruit_set_success/(fruit.Fruit_set_success+fruit.Fruit_set_failure)
    seeds = pd.read_csv(folder/"raspseed.csv", sep=";").query("Plant_sp == 'RB'").groupby("Treatment").Seed_set.mean()
    rows = [
        ("Rubus idaeus", "SmallFruit_2025", "10.1111/jen.70019", "fruit set", "Figshare 29930195, successes/all flowers", fruit_rate.Day,fruit_rate["Day+Night"],True),
        ("Rubus idaeus", "SmallFruit_2025", "10.1111/jen.70019", "seed set", "Figshare 29930195, seeds per fruit", seeds.Day,seeds["Day+Night"],False),
        ("Fagopyrum esculentum", "Nagano_2024", "10.1007/s11829-024-10117-x", "seed set", "Fig 5 medians, y=521 for 0 and y=27 for 0.6; bar centres y=294 and 223", (521-294)*0.6/494,(521-223)*0.6/494,True),
        ("Carica papaya", "Otieno_2026", "10.64898/2026.07.17.739175", "fruit set", "Preprint Table S1 adjusted fruit count ratio, open/day=1.58; normalized to open=1",1/1.58,1.,True),
    ]
    return pd.DataFrame(rows, columns=["plant","study_ID","source",METRIC,"basis","day_value","open_value","used"])


def analyze(inputs, signed, folder, out):
    out.mkdir(parents=True, exist_ok=True)
    cohort = set(pd.read_csv(ROOT/"results/moth_pollen_evidence/shared_544_evidence.csv").plant)
    assert len(cohort) == 544
    pairs, screened = paired_experiments(folder,inputs,cohort)
    pairs = pd.concat([pairs,recent(folder)],ignore_index=True)
    assert pairs.open_value.gt(0).all() and pairs.day_value.ge(0).all()
    pairs["loss"] = 1-pairs.day_value/pairs.open_value
    pairs.to_csv(out/"experiments.csv",index=False)
    per_study = pairs[pairs.used].groupby(["plant","study_ID"]).loss.mean()
    traits = per_study.groupby("plant").agg(["mean","count"]).rename(columns={"mean":"loss","count":"studies"})
    assert set(traits.index) <= cohort
    screen = pd.DataFrame(index=pd.Index(sorted(cohort),name="plant"))
    screen["indexed_primary_studies"] = screened.groupby("plant").study_ID.agg(lambda s: "|".join(sorted(set(s))))
    screen["estimated"] = screen.index.isin(traits.index)
    screen.to_csv(out/"coverage_544.csv")
    sites, digest, hosts, _ = primary(inputs,signed)
    pos = hosts[hosts.wcvp_accepted_name.isin(traits.index)][["key","wcvp_accepted_name"]].drop_duplicates()
    present = pd.crosstab(pos.key,pos.wcvp_accepted_name).reindex(index=sites.key,columns=traits.index,fill_value=0)
    assert present.max().max()==1 and len(sites)==10698 and sites.region5.nunique()==78
    def estimate(v):
        return fit(sites.assign(exposure_moth=v),"exposure","moth",balance=True)
    for plant in traits.index:
        hits = present[plant].to_numpy()
        model = estimate(hits)
        traits.loc[plant,"gap"] = model[GAP]
        for land in ["island","mainland"]:
            traits.loc[plant,land+"_positive_site_years"] = hits[sites.land_class.eq(land)].sum()
    centered = traits.loss-traits.loss.mean()
    weights = centered/centered.pow(2).sum()
    result = estimate(present.to_numpy() @ weights.to_numpy())
    np.testing.assert_allclose(result[GAP],np.polyfit(traits.loss,traits.gap,1)[0],atol=1e-12)
    traits.to_csv(out/"species.csv")
    positive = sites.loc[present.any(axis=1).to_numpy()]
    result.update(species=len(traits), studies=len(per_study),target_keys_sha256=digest,
                  regions_with_species=positive.region5.nunique(),
                  regions_positive_both_sides=int(positive.groupby("region5").land_class.nunique().eq(2).sum()),
                  positive_site_years=positive.land_class.value_counts().to_dict())
    (out/"contrast.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["inputs","signed","sources","out"]:
        parser.add_argument(arg,type=Path)
    args = parser.parse_args()
    analyze(args.inputs,args.signed,args.sources,args.out)
