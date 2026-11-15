#!/usr/bin/env python3
"""Li 2026 original 79 plants, shared site x species x flowering-round support.

NOT a new independent experiment, not factorial geometry-only causal effect.
Data SHA pin: author Git blob b616aca5beeba1457549841477d61ca4f796a147.
A stratum must have authentic L and C rows with measured seeds and Starting>0.
"""
from __future__ import annotations
import hashlib, json, itertools, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from io import BytesIO
ROOT=Path(__file__).resolve().parents[2]
URL="https://raw.githubusercontent.com/Dongboli/experimental-data/main/Artificial%20lights%20at%20night/seedset.csv"
BLOB="b616aca5beeba1457549841477d61ca4f796a147"
OUT=ROOT/"results/alan_Li2026_matched_site_species_round_seed_bud_20261115.json"
STRATA=("replicate","species","round")
def fetch():
    req=urllib.request.Request(URL,headers={"User-Agent":"ALAN-Li2026-authentic-stratum-gate"})
    with urllib.request.urlopen(req,timeout=50) as r:b=r.read(50000)
    sha=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()
    if sha!=BLOB:raise ValueError("Author Git object identity mismatch")
    return pd.read_csv(BytesIO(b))
def prepare(raw):
    if len(raw)!=79: raise AssertionError("Source observation count drift")
    need={"Treatment","replicate","species","round","Starting","Fruit","Seeds","pollinated"}
    if not need<=set(raw.columns):raise AssertionError("Missing author columns")
    x=raw.copy()
    x["Starting"]=pd.to_numeric(x.Starting,errors="coerce")
    x["Fruit"]=pd.to_numeric(x.Fruit,errors="coerce")
    x["Seeds"]=pd.to_numeric(x.Seeds,errors="coerce")
    x["replicate"]=pd.to_numeric(x.replicate,errors="coerce")
    x["round"]=pd.to_numeric(x["round"],errors="coerce")
    if set(x.Treatment.dropna())!={"L","C","P"}:raise AssertionError("Treatment arms drift")
    if set(x.replicate.dropna())!={1,2,3}:raise AssertionError("Three original fields drift")
    if set(x.species.dropna())!={"ML","SN","ZC"}:raise AssertionError("Species drift")
    if (x.Seeds.notna().sum()!=78):raise AssertionError("Missingness drift")
    if ((x.Starting.dropna()<0).any() or (x.Fruit.dropna()<0).any() or (x.Seeds.dropna()<0).any()):
        raise AssertionError("Negative raw source plant outcomes")
    return x
def matched(raw):
    x=prepare(raw)
    # Correct exclusion does not mark unmeasured seeds zero and does not
    # assume blooming plants with Starting<=0 had any true flower opportunity.
    elig=x[x.Seeds.notna() & x.Fruit.notna() & x.Starting.gt(0) & x.replicate.notna()
         & x["round"].notna()].copy()
    pairs=[]
    for key,g in elig.groupby(list(STRATA),sort=True,dropna=False):
        if not {"C","L"}<=set(g.Treatment):continue
        arms={}
        for light in ("C","L"):
            a=g[g.Treatment==light]
            arms[light]={"n_plants":int(len(a)),"buds":float(a.Starting.sum()),
                         "fruits":float(a.Fruit.sum()),"seeds":float(a.Seeds.sum()),
                         "seed_positive":int((a.Seeds>0).sum())}
        if min(arms["C"]["buds"],arms["L"]["buds"])<=0:continue
        dark=arms["C"];lit=arms["L"]
        rec={"field":int(key[0]),"species":str(key[1]),"round":int(key[2]),
             "dark_n_plants":dark["n_plants"],"lit_n_plants":lit["n_plants"],
             "dark_buds":dark["buds"],"lit_buds":lit["buds"],
             "dark_seeds":dark["seeds"],"lit_seeds":lit["seeds"],
             "dark_fruits":dark["fruits"],"lit_fruits":lit["fruits"],
             "dark_seed_per_bud":dark["seeds"]/dark["buds"],
             "lit_seed_per_bud":lit["seeds"]/lit["buds"],
             "dark_fruit_per_bud":dark["fruits"]/dark["buds"],
             "lit_fruit_per_bud":lit["fruits"]/lit["buds"],
             "dark_seed_positive":dark["seed_positive"],"lit_seed_positive":lit["seed_positive"]}
        rec["seed_per_bud_diff"]=rec["lit_seed_per_bud"]-rec["dark_seed_per_bud"]
        rec["fruit_per_bud_diff"]=rec["lit_fruit_per_bud"]-rec["dark_fruit_per_bud"]
        pairs.append(rec)
    return pairs,{"source_missing_seeds_by_original_treatment":
                    {str(k):int(v) for k,v in x[x.Seeds.isna()].Treatment.value_counts().items()},
                  "matched_strata_dark_seed_zero":
                    sum(p["dark_seeds"]==0 for p in pairs),
                  "matched_strata_n_ge_2_plants_both_arms":
                    sum(p["dark_n_plants"]>=2 and p["lit_n_plants"]>=2 for p in pairs),
                  "matched_strata_one_plant_either_arm":
                    sum(p["dark_n_plants"]==1 or p["lit_n_plants"]==1 for p in pairs),
                  "author_observations":len(x),"known_seed_observations":int(x.Seeds.notna().sum()),
                  "eligible_with_start_buds":len(elig),"source_fields":3,"matched_strata":len(pairs),
                  "unmatched_or_missing_not_zero":True}
def sign_p(nplus,nminus):
    from math import comb
    n=nplus+nminus
    if not n:return None
    return float(min(1.,2*sum(comb(n,i) for i in range(min(nplus,nminus)+1))/2**n))
def stats(pairs):
    if not pairs:return {"status":"no_matched_strata"}
    out={}
    for measure in ("seed_per_bud_diff","fruit_per_bud_diff"):
        vals=[float(p[measure]) for p in pairs]
        byfield={str(k):float(np.mean([p[measure] for p in pairs if p["field"]==k]))
                 for k in [1,2,3] if any(p["field"]==k for p in pairs)}
        q=np.array(list(byfield.values()),float)
        signs=(int(sum(z>0 for z in vals)),int(sum(z<0 for z in vals)),int(sum(z==0 for z in vals)))
        field_signs=(int((q>0).sum()),int((q<0).sum()),int((q==0).sum()))
        out[measure]={"equal_stratum_mean_difference":float(np.mean(vals)),
            "strata_positive_negative_zero":signs,
            "true_field_equal_mean_difference":float(np.mean(q)),
            "field_equal_differences":byfield,
            "field_signs_positive_negative_zero":field_signs,
            "exact_3_field_two_sided_sign_probability_ONLY_if_exchangeable":
                sign_p(field_signs[0],field_signs[1]),
            "not_a_confirmed_randomization_pvalue":True}
    by_species={}
    for s in sorted(set(p["species"] for p in pairs)):
        group=[p for p in pairs if p["species"]==s]
        by_species[s]={"matched_strata":len(group),
                       "equal_stratum_mean_seed_per_bud_difference":
                        float(np.mean([p["seed_per_bud_diff"] for p in group]))}
    out["by_species"]=by_species
    for drop in sorted(by_species):
        left=[p for p in pairs if p["species"]!=drop]
        out["leave_one_species_out_"+drop]={"n_strata":len(left),
            "equal_stratum_mean_seed_per_bud_difference":
               float(np.mean([p["seed_per_bud_diff"] for p in left])) if left else None}
    # Reweight available matched strata: aggregate within matched source,
    # rather than all source rows which may have different taxon/round support.
    for stage in ("seeds","fruits"):
        num=sum(p["lit_"+stage] for p in pairs)/sum(p["lit_buds"] for p in pairs)
        den=sum(p["dark_"+stage] for p in pairs)/sum(p["dark_buds"] for p in pairs)
        out[stage+"_matched_strata_pooled_rate_ratio"]=float(num/den) if den>0 else None
    out["not_a_model_based_ALAN_meta_effect"]=True
    return out
def main():
    pairs,audit=matched(fetch())
    supported2=[p for p in pairs if p["dark_n_plants"]>=2 and p["lit_n_plants"]>=2]
    out={"schema":"ALAN_Li2026_true_field_x_species_x_round_common_support_v2",
         "original_author_git_blob":BLOB,"original_doi":"10.64898/2026.07.10.737511",
         "source_audit":audit,"matched_strata":pairs,
         "stats":stats(pairs),
         "minimum_two_plants_per_arm_matched_strata":
              {"n":len(supported2),"results":stats(supported2),
               "descriptive_sensitivity_NOT_independent_field_replications":True},
         "limits":["Only three independently manipulated fields; plant rows are not experimental light replicates",
                   "Equal-stratum difference is within-field descriptive and is not a field-randomization estimate",
                   "Round/species support is incomplete and can restrict the sample; unmatched plant rows excluded",
                   "All treatments share field/plant reuse concerns; site x round strata are not independent replicates",
                   "Initial buds are opportunity proxies; not original fully exposed mature flowers",
                   "No source illumination geometry-only causal effect or moth pollen-flow mediation test",
                   "Not an ALAN cross-study meta-analysis, not island-dependent study"]}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("ALAN_LI_REAL_MATCHED="+json.dumps({k:v for k,v in out.items() if k!="matched_strata"},ensure_ascii=False))
if __name__=="__main__":main()
