"""Regional macroecological trait contrasts and held-out-region prediction."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import brier_score_loss

from alan_globi_v011_same_fiveyear import taxonomy


def flag(series, states):
    return series.isin(states).where(series.notna()).astype(float)


def run(work):
    root = Path(__file__).resolve().parents[2] / "results"
    src = root / "regional_light_traits"
    out = src / "mechanisms"
    out.mkdir(exist_ok=True)
    traits = pd.read_csv(src / "species_traits.csv").set_index("sp")
    name = traits.accepted_name
    direct = pd.read_csv(work / "direct_trait_replay/direct_species_trait_ledger.csv.gz", low_memory=False)
    direct = direct.loc[direct.quality.isin(["medium", "high"]) & direct.resolution_status.eq("resolved")
                        & direct.multistate.eq(False) & direct.normalized_value.notna()
                        & ~direct.normalized_value.str.contains("|", regex=False, na=False)]
    axes = direct.pivot(index="accepted_species", columns="trait_name", values="normalized_value")
    a = axes.reindex(name).set_axis(traits.index)
    t = pd.DataFrame(index=traits.index)
    t["annual"] = traits.lifeform.str.contains(r"\bannual\b", case=False).astype(float)
    t["white_tubular"] = (flag(a.flower_primary_color, ["white"]) *
                          flag(a.floral_form, ["tubular", "salverform", "funnel_trumpet"]))
    t["white_tubular_deep"] = t.white_tubular * flag(a.tube_depth_class, ["deep"])
    t["self_incompatible"] = flag(a.self_incompatibility.where(a.self_incompatibility.isin(["SC","SI"])), ["SI"])
    t["autonomous_selfing"] = flag(a.autonomous_selfing_capacity, ["autonomous"])
    t["predominantly_selfing"] = flag(a.mating_system, ["predominantly_selfing", "obligate_selfing"])
    mating = pd.read_csv(work / "mechanism_originals/Database_BioGeo_mating.csv", encoding="latin1")
    pollen = pd.read_csv(work / "pollen/gl.csv", encoding="latin1")
    mating["raw"] = mating.genus.str.strip() + " " + mating.species.str.strip()
    pollen["raw"] = pollen.Species_accepted_names.str.replace("_", " ", regex=False).str.strip()
    cross = taxonomy(work / "v011", set(pd.concat([mating.raw,pollen.raw]).dropna()), set())
    matched = cross.loc[cross.resolution.eq("resolved"),"accepted_name"]
    mating["selfing"] = (1-mating["mean.tm"]).where(mating["mean.tm"].between(0,1))
    t["genetic_selfing_fraction"] = name.map(mating.groupby(mating.raw.map(matched)).selfing.median())
    t["pollen_limitation"] = name.map(pollen.groupby(pollen.raw.map(matched)).PL_Effect_Size.median())
    reviewed = traits.transport_review_status.copy()
    updates = pd.read_csv(root / "mechanism_exemplars/transport_updates.csv")
    reviewed.loc[name.isin(updates.plant)] = "documented"
    t["moth_pollen_documented"] = flag(reviewed.where(reviewed.isin(["documented","unconfirmed"])), ["documented"])
    t["night_exclusion_loss"] = traits.night_exclusion_loss
    flags = pd.read_csv(next((work / "v011").glob("*flags/*csv.gz")))
    links = pd.read_csv(work / "v011/moth_name_crosswalk.csv")
    positive = flags.loc[flags.non_inat_visits_moth_species_insect.eq(1),"plant_name"]
    documented = links.loc[links.resolution.eq("resolved") & links.name.isin(positive),"accepted_name"]
    t["moth_visit_documented"] = flag(name, documented)
    data = pd.read_csv(src / "all_regions.csv.gz").query("blocks >= 3").join(t,on="sp")
    data["band"] = pd.cut(data.blocks,[2,5,10,25,np.inf], labels=False)
    data["negative"] = data.label.eq("neg")
    assert data.loc[data.negative,"sp"].nunique() == 696
    predator_chains = pd.read_csv(root / "species_interpretation/land_negative_predator_chains.csv")
    hits = data.loc[data.negative, ["sp","region","land_class"]].join(name,on="sp").merge(
        predator_chains,left_on="accepted_name",right_on="plant")
    hits[["sp","region","land_class","plant","moth","predator"]].drop_duplicates().to_csv(
        out/"candidate_plant_chain_name_matches.csv",index=False)
    rows = []
    for land, frame in [("both",data), *data.groupby("land_class")]:
        for col in t:
            candidates = frame.loc[frame.negative]
            known = candidates.dropna(subset=[col])
            ref = frame.loc[~frame.negative].groupby(["region","band"])[col].mean()
            pairs = known.join(ref.rename("ref"),on=["region","band"]).dropna(subset=["ref"])
            diff = pairs[col]-pairs.ref
            loo = [(diff[pairs.region.ne(r)]).mean() for r in pairs.region.unique()]
            rows.append(dict(land=land,trait=col,total=candidates.sp.nunique(),known=known.sp.nunique(),
                             paired_species=pairs.sp.nunique(),pairs=len(pairs),regions=pairs.region.nunique(),
                             candidate=pairs[col].mean(),reference=pairs.ref.mean(),difference=diff.mean(),
                             loo_min=pd.Series(loo,dtype=float).min(),loo_max=pd.Series(loo,dtype=float).max()))
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "matched_traits.csv",index=False)
    data = data.join(traits.wcvp_family, on="sp")
    strata = ["region","band","annual","wcvp_family"]
    adjusted = []
    for land, frame in data.groupby("land_class"):
        for col in ["autonomous_selfing","predominantly_selfing","self_incompatible"]:
            ref = frame.loc[~frame.negative].groupby(strata)[col].mean().rename("ref")
            pairs = frame.loc[frame.negative].join(ref,on=strata).dropna(subset=[col,"ref"])
            adjusted.append(dict(land=land,trait=col,species=pairs.sp.nunique(),pairs=len(pairs),
                                 regions=pairs.region.nunique(),candidate=pairs[col].mean(),
                                 reference=pairs.ref.mean(),difference=(pairs[col]-pairs.ref).mean()))
    pd.DataFrame(adjusted).to_csv(out/"family_annual_matched.csv",index=False)
    # Compare prediction for new regions on the same species-region panel and folds.
    # Fixed unknown-state coding; missingness can itself predict recording practices.
    x = pd.DataFrame({"log_blocks":np.log(data.blocks),"island":data.land_class.eq("island").astype(int)})
    folds = []
    for r, group in data.groupby("region"):
        test = data.region.eq(r)
        row = dict(region=r,n=len(group),negative=int(group.negative.sum()))
        for model,features in [("support",x),("with_traits",x.join(data[t.columns].fillna(-999)))]:
            rf = RandomForestClassifier(n_estimators=100,min_samples_leaf=20,max_features=1.0,random_state=42,n_jobs=3)
            rf.fit(features[~test],data.negative[~test])
            pred = rf.predict_proba(features[test])[:,1]
            row[model] = brier_score_loss(group.negative,pred)
        folds.append(row)
    folds = pd.DataFrame(folds)
    folds.to_csv(out / "held_out_regions.csv",index=False)
    fig, ax = plt.subplots(figsize=(8,4.6))
    binary = ["annual","white_tubular","self_incompatible","autonomous_selfing","predominantly_selfing","moth_pollen_documented","moth_visit_documented"]
    for land,offset,color in [("island",-.13,"#256c79"),("mainland",.13,"#a55b38")]:
        z=summary[summary.land.eq(land)].set_index("trait").reindex(binary)
        ax.scatter(100*z.difference,np.arange(len(z))+offset,label=land,color=color)
    ax.set_yticks(range(len(binary)),["Annual life form","White, tubular flower","Self-incompatible","Autonomous selfing capacity","Predominantly/obligately selfing","Moth pollen transport documented*","Moth visit documented"])
    ax.axvline(0,color="0.65",lw=1); ax.invert_yaxis(); ax.legend()
    ax.set_xlabel("Candidate minus matched reference (percentage points)")
    fig.text(.02,.01,"* Among reviewed species. Coverage and region counts vary by trait; see matched_traits.csv.",fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1)); fig.savefig(out/"traits_by_land.png",dpi=180); plt.close(fig)


if __name__ == "__main__":
    run(Path(sys.argv[1]))
