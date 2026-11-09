#!/usr/bin/env python3
"""EXACT global 93,209,351 original-photo moth-positive opportunities by iNat taxon.

Use FCP's ORIGINAL signed 42,111-species 2026-09-27 source count CSV; each
inat_taxon_id has the true number of unique original photo IDs in the full
93m archive (NOT 93m actual photographs with flower petals/ROI9).
Join WCVP exact accepted IDs / names and signed published GloBI study-linked
flower interactions at the accepted binomial. GloBI moth-positive is a
VISIT-RECORD criterion; no record is UNKNOWN, not moth absence. A GloBI
'pollinates' relationship is not proof of nocturnal effective pollen deposition.
The Walton primary study supports pollen carriage by night-sampled moths,
not all-site moth dependence. Never infer native species, local visitation,
actual flower blooming, or original RGB availability from per-species totals.

This source-based exact count is far faster and more faithful than repeating
all 224 original 93m tier archives just to count photo IDs; later use tier
archives ONLY to recover individual IDs, years, geography, license and bytes.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import pandas as pd

SCHEMA="ALAN_full93209351_original_photo_WCVP_GloBI_moth_evidence_coverage_v1"
PINS={
 "counts":"ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9",
 "wcvp":"dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2",
 "accepted":"9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3",
 "globi":"a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994",
}
N_PHOTO=93209351
N_SOURCE=42111
N_GLOBI=61967
EVIDENCE={
 "moth_positive":"moth_visit_positive",
 "moth_pollinates_db":"moth_pollinates_db",
 "walton2020_night_pollen_carriage":"walton2020_night_pollen_carriage",
 "butterfly_positive":"butterfly_positive",
}

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_GLOBAL_93M_SOURCE_MOTH: "+msg)

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for q in iter(lambda:f.read(4<<20),b""):h.update(q)
    return h.hexdigest()

def evidence(g):
    require({"plant_name","guild","relation","study_title"}<=set(g),
            "Missing GloBI study-linked moth source fields")
    require(set(g.guild)<=set(("moth_family_visit","butterfly_family_visit")) and
            set(g.relation)<=set(("visitsFlowersOf","pollinates")),
            "Wrong orientation or non-flower GloBI interaction")
    t=g.assign(moth_visit_positive=g.guild.eq("moth_family_visit"),
         moth_pollinates_db=(g.guild.eq("moth_family_visit")&
                             g.relation.eq("pollinates")),
         walton2020_night_pollen_carriage=(g.guild.eq("moth_family_visit")&
                             g.study_title.str.contains("rsbl.2019.0877",
                                 case=False,regex=False,na=False)),
         butterfly_positive=g.guild.eq("butterfly_family_visit"))
    return t.groupby("plant_name",observed=True)[list(EVIDENCE.values())].any().reset_index()

def compute_counts(counts,wcvp,accepted,globi,*,strict_counts=True):
    assert counts.inat_taxon_id.is_unique
    n=counts.copy()
    n["n_unique_source_photo_ids"]=pd.to_numeric(
        n.n_unique_source_photo_ids,errors="raise").astype("int64")
    require(n.n_unique_source_photo_ids.ge(0).all(),
            "Negative original photo IDs for an iNaturalist species")
    if strict_counts:
        require(len(n)==N_SOURCE and n.n_unique_source_photo_ids.sum()==N_PHOTO,
                "Original 42,111 taxon / 93,209,351 photo ID source counts drift")
    w=wcvp[["inat_taxon_id","wcvp_match_state",
            "wcvp_accepted_plant_name_id"]].copy()
    require(w.inat_taxon_id.is_unique,"Original taxa source not unique")
    w=w.loc[w.wcvp_match_state.eq("unique_exact_name_to_accepted_id")].copy()
    a=accepted[["wcvp_accepted_plant_name_id","wcvp_accepted_name"]].drop_duplicates()
    require(a.wcvp_accepted_plant_name_id.is_unique,
            "Accepted WCVP botanical name inconsistent across countries")
    x=n.merge(w,on="inat_taxon_id",how="left",validate="one_to_one")
    x=x.merge(a,on="wcvp_accepted_plant_name_id",how="left",validate="many_to_one")
    f=evidence(globi)
    x=x.merge(f,left_on="wcvp_accepted_name",right_on="plant_name",
              how="left",validate="many_to_one")
    require(len(x)==len(n) and x.inat_taxon_id.is_unique and
            n.n_unique_source_photo_ids.sum()==x.n_unique_source_photo_ids.sum(),
            "Evidence join multiplied or lost any original photo IDs")
    for key in EVIDENCE.values():
        x[key]=x[key].fillna(False).astype(bool)
    not_matched=x.wcvp_accepted_name.isna()
    out={
      "schema":SCHEMA,
      "source_photo_count":int(n.n_unique_source_photo_ids.sum()),
      "source_distinct_inat_taxa":len(n),
      "source_distinct_valid_WCVP_accepted_species":int(
          x.wcvp_accepted_name.nunique()),
      "photo_IDs_in_source_taxa_without_verified_WCVP_accepted_name":
          int(x.loc[not_matched,"n_unique_source_photo_ids"].sum()),
      "source_taxa_not_verified_WCVP_accepted_name":int(not_matched.sum()),
      "source_GloBI_study_linked_interactions":len(globi),
      "GloBI_positive_moth_names_in_source":int(
          f.moth_visit_positive.sum()),
      "evidence_tiers":{},
      "positive_both_moth_and_butterfly_original_photo_IDs":int(
          x.loc[x.moth_visit_positive & x.butterfly_positive,
                "n_unique_source_photo_ids"].sum()),
      "moth_positive_without_recorded_butterfly_original_photo_IDs":int(
          x.loc[x.moth_visit_positive & ~x.butterfly_positive,
                "n_unique_source_photo_ids"].sum()),
      "both_guild_source_positive_is_NOT_direct_local_moth_or_butterfly_dominance":True,
      "photographic_source_iNaturalist_study_record_ascertainment_risk":True,
      "unknown_GloBI_record_cannot_be_called_moth_absence":True,
      "source_exhausts_true_photo_ID_OPPORTUNITIES_not_true_flower_ROI_pixels":True,
      "not_all_wild_flora_or_native_botanical_colour":True,
    }
    source_by_taxon=[]
    for label,key in EVIDENCE.items():
        d=x.loc[x[key]].copy()
        out["evidence_tiers"][label]={
            "original_distinct_photo_IDs":int(d.n_unique_source_photo_ids.sum()),
            "unique_original_iNat_taxa":int(d.inat_taxon_id.nunique()),
            "unique_WCVP_accepted_species":int(d.wcvp_accepted_name.nunique()),
            "fraction_of_original_93m_source_photo_IDs":float(
                d.n_unique_source_photo_ids.sum()/
                max(1,n.n_unique_source_photo_ids.sum())),
            "512photo_worker_groups_needed_to_measure_all_original_RGB":int(
                math.ceil(d.n_unique_source_photo_ids.sum()/512)),
            "disjoint_from_other_evidence_tiers":False,
        }
        ranked=d.sort_values("n_unique_source_photo_ids",ascending=False)
        out["evidence_tiers"][label]["ten_most_photographed_source_taxa"]=[
             {"accepted_name":str(t.wcvp_accepted_name),
              "original_photo_IDs":int(t.n_unique_source_photo_ids)}
             for t in ranked.head(10).itertuples(index=False)]
        # Publish only aggregate source counts, original species name/accepted
        # ID; omit original individual photo/observation/observer coordinates.
        source_by_taxon.append(d[[
            "inat_taxon_id","species","wcvp_accepted_plant_name_id",
            "wcvp_accepted_name","n_unique_source_photo_ids"]].assign(evidence_tier=label))
    return out,pd.concat(source_by_taxon,ignore_index=True)

def run(paths,out):
    out=Path(out);require(not out.exists(),"Refuse to overwrite source-verified 93m moth counts")
    for k,p in paths.items():
        require(sha(p)==PINS[k],"Frozen original source checksum changed: "+k)
    counts=pd.read_csv(paths["counts"],dtype=str,keep_default_na=False)
    wcvp=pd.read_csv(paths["wcvp"],compression="gzip",dtype=str,keep_default_na=False)
    acc=pd.read_csv(paths["accepted"],compression="gzip",dtype=str,keep_default_na=False,
        usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    globi=pd.read_csv(paths["globi"],compression="gzip",dtype=str,keep_default_na=False)
    require(len(globi)==N_GLOBI,"GloBI original 61,967 study links not source-complete")
    summary,detail=compute_counts(counts,wcvp,acc,globi)
    require(summary["evidence_tiers"]["moth_positive"]["original_distinct_photo_IDs"]
             ==34346428 and
            summary["evidence_tiers"]["moth_pollinates_db"]["original_distinct_photo_IDs"]
             ==13262462 and
            summary["evidence_tiers"]["walton2020_night_pollen_carriage"]["original_distinct_photo_IDs"]
             ==3738985,
            "Independent old source counts/accepted name MOTH evidence changed")
    summary["source_SHA256"]=PINS
    out.mkdir(parents=True)
    detail_file=out/"accepted_species_original93m_photo_counts_by_moth_evidence.csv"
    detail.to_csv(detail_file,index=False)
    summary["result_species_tier_counts_sha256"]=sha(detail_file)
    (out/"ALAN_full93m_original_photo_moth_positive_opportunity_count.json").write_text(
        json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print("ALAN_SOURCE_AUTHENTICATED_FULL93M_MOTH_RESTRICTED_PHOTO_ID_COUNTS",
          json.dumps({k:v for k,v in summary.items() if k in (
           "source_photo_count","source_distinct_inat_taxa",
           "photo_IDs_in_source_taxa_without_verified_WCVP_accepted_name",
           "positive_both_moth_and_butterfly_original_photo_IDs",
           "moth_positive_without_recorded_butterfly_original_photo_IDs",
           "evidence_tiers")},sort_keys=True),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in PINS:p.add_argument("--"+k,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run({k:getattr(args,k) for k in PINS},args.out)
if __name__=="__main__":
    main()
