#!/usr/bin/env python3
"""Original Li 2026 triplet L/P/C within same field x species x flowering round.

This is a descriptive precision-source audit, not a causal geometry-only test.
The original 79 plant rows are verified by author Git blob SHA in imported fetch().
"""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import json
import numpy as np
from scripts.analysis.alan_Li2026_matched_field_taxon_round_20261115 import (
    ROOT, BLOB, fetch, prepare, STRATA
)
OUT=ROOT/"results/alan_Li2026_original_triple_arm_light_geometry_and_dark_floor_20261115.json"

def groups(raw):
    x=prepare(raw)
    elig=x.loc[x.Seeds.notna() & x.Fruit.notna() & x.Starting.gt(0)
               & x.replicate.notna() & x["round"].notna()].copy()
    q=[]
    for k,g in elig.groupby(list(STRATA),sort=True,dropna=False):
        if not {"L","P","C"} <= set(g.Treatment):continue
        arms={}
        for light in ("L","P","C"):
            a=g[g.Treatment==light]
            arms[light]={"plants":int(len(a)),"buds":float(a.Starting.sum()),
                         "fruits":float(a.Fruit.sum()),"seeds":float(a.Seeds.sum()),
                         "seed_positive":int((a.Seeds>0).sum())}
        if min(v["buds"] for v in arms.values())<=0:continue
        row={"field":int(k[0]),"species":str(k[1]),"round":int(k[2]),
             "arms":arms,"dark_zero":arms["C"]["seeds"]==0,
             "point_zero":arms["P"]["seeds"]==0}
        for out in ("seeds","fruits"):
            for t in ("L","P","C"):
                row[out+"_"+t]=arms[t][out]/arms[t]["buds"]
            for t,u in (("L","C"),("P","C"),("L","P")):
                row[out+"_"+t+"minus"+u]=row[out+"_"+t]-row[out+"_"+u]
        q.append(row)
    return q,{"original_plants":len(x),"original_seeds_observed":int(x.Seeds.notna().sum()),
              "source_fields":int(x.replicate.nunique()),"three_arm_strata":len(q),
              "unmatched_strata_do_not_qualify":True}
def signs(vals):
    return {"positive":sum(v>0 for v in vals),
            "negative":sum(v<0 for v in vals),
            "zero":sum(v==0 for v in vals)}
def summary(q):
    if not q:return {"n":0,"not_estimated":True}
    out={"n":len(q),"n_dark_zero":sum(z["dark_zero"] for z in q),
         "n_point_zero":sum(z["point_zero"] for z in q),
         "fields":len(set(z["field"] for z in q)),
         "species":sorted(set(z["species"] for z in q)),
         "field_signs_and_differences":{}}
    for outcome in ("seeds","fruits"):
        for a,b in (("L","C"),("P","C"),("L","P")):
            label=outcome+"_"+a+"minus"+b
            vals=[z[label] for z in q]
            site={str(f):float(np.mean([z[label] for z in q if z["field"]==f]))
                  for f in sorted(set(z["field"] for z in q))}
            count=signs(vals)
            out[label]={"strata_signs":count,
                         "equal_strata_mean_difference":float(np.mean(vals)),
                         "field_equal_mean_of_stratum_differences":
                             float(np.mean(list(site.values()))),
                         "field_means":site,
                         "field_signs":signs(list(site.values())),
                         "not_independent_strata":True}
        for a,b in (("L","C"),("P","C"),("L","P")):
            denom=sum(z["arms"][b]["buds"] for z in q)
            num=sum(z["arms"][a]["buds"] for z in q)
            r_b=sum(z["arms"][b][outcome] for z in q)/denom
            r_a=sum(z["arms"][a][outcome] for z in q)/num
            out[outcome+"_"+a+"over"+b+"_pooled_rate_ratio"] = float(r_a/r_b) if r_b>0 else None
    out["no_meta_SE_or_randomization_p"]=True
    return out
def analyze(raw):
    q,orig=groups(raw)
    common2=[z for z in q if all(z["arms"][t]["plants"]>=2 for t in ("L","P","C"))]
    dark_positive=[z for z in q if not z["dark_zero"]]
    return {
        "schema":"ALAN_Li2026_original_triarm_matched_zero_floor_competitor_v1",
        "source_doi":"10.64898/2026.07.10.737511",
        "author_git_blob":BLOB,
        "source":orig,"all_three_arm_matched":summary(q),
        "all_three_arm_matched_source_strata":q,
        "dark_nonzero_only_POSTHOC_SELECTION":summary(dark_positive),
        "at_least_two_plants_per_each_L_P_C_arm":summary(common2),
        "interpretation":{
            "L":"25m linear array, approx 26 lamps",
            "P":"isolated point lamps, approx 4 lamps",
            "C":"unlit",
            "L_vs_P_is_not_equal_lux_equal_area_geometry_only":True,
            "conditioning_on_dark_zero_is_outcome_selected":True,
            "three_fields_only_not_11_independent_lamp_replicates":True,
            "no_moth_or_stigma_pollen_data":True,
            "no_cross_study_ALAN_meta_beta":True
        }
    }
def main():
    data=analyze(fetch())
    assert data["source"]["original_plants"]==79
    assert data["source"]["original_seeds_observed"]==78
    assert data["source"]["source_fields"]==3
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("ALAN_LI_TRIARM_SOURCE="+json.dumps({k:v for k,v in data.items()
          if k!="all_three_arm_matched_source_strata"},ensure_ascii=False))
if __name__=="__main__":main()
