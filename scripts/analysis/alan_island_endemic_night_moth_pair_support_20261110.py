#!/usr/bin/env python3
"""First empirical feasibility test of island-endemic nocturnal-pollination pairs.

Original data: signed ALAN 93,209,351 photograph-ID source as 16 group SQLite
WCVP-accepted taxon/island/site/year tables; real PCNL 1992-2024 (prior only);
Kendall & Nicholson (2025) original day/night exclusion study/trait source.

A documented NIGHT-blooming plant is NOT proof of moth dependence.
A hawkmoth-effective island plant is NOT evidence that a sister lineage has
undergone a documented day-to-night evolutionary switch.
A species observed on an island is NOT automatically endemic or wild there.
This program is only a provenance-checked opportunity census, never a 3-way
effect fit or GloBI biological-absence imputation.
"""
from __future__ import annotations
import argparse,hashlib,json,itertools,sqlite3
from pathlib import Path
import numpy as np
import pandas as pd

SITE_SHA="440a8313d71168e0a4dfd2f41a7bef45623abb068b07e3282026f1c35dccbd2c"
STUDY_SOURCE_BLOB="18381e550a7e405fe4f198909f0bb381bbecf1e0"
KEY=["island","latbin","lonbin","year"]
YEARS={2019,2021,2023,2024,2025}
SCHEMA="ALAN_island_literature_nocturnal_effective_moth_day_control_pair_support_v1"
# Seed set only. Not an exhaustive global island endemic catalogue.
SEEDS={
 "Ipomoea habeliana":dict(archipelago="Galapagos",evidence="effective_hawkmoth_documented; facultative_selfing",source="doi:10.1111/j.1095-8339.2009.00963.x",night_moth_effective=True,evolutionary_day_to_night_shift_proven=False,endemic_range="confirmed_Galapagos_archipelago"),
 "Brighamia insignis":dict(archipelago="Hawaii",evidence="moth_syndrome_only; no_verified_effective_moth; extinct_in_wild",source="doi:10.1016/j.flora.2019.151448",night_moth_effective=False,evolutionary_day_to_night_shift_proven=False,endemic_range="confirmed_Hawaii_archipelago"),
 "Digitalis canariensis":dict(archipelago="Canary_Islands",evidence="effective_passerine_bird; accepted_name_for_Isoplexis_canariensis",source="doi:10.3732/ajb.0800055",night_moth_effective=False,evolutionary_day_to_night_shift_proven=False,endemic_range="confirmed_Canary_archipelago"),
 "Canarina canariensis":dict(archipelago="Canary_Islands",evidence="effective_passerine_bird; floral_herbivory_possible",source="PMID:18931830; PMID:21900611",night_moth_effective=False,evolutionary_day_to_night_shift_proven=False,endemic_range="confirmed_Canary_archipelago"),
 "Echium simplex":dict(archipelago="Tenerife",evidence="flying_insects_effective; not_proven_strict_day_only",source="doi:10.1111/plb.12926",night_moth_effective=False,evolutionary_day_to_night_shift_proven=False,endemic_range="confirmed_Tenerife"),
 "Cordia lutea":dict(archipelago="Galapagos_and_mainland",evidence="native_but_not_endemic_Galapagos; diurnal_and_nocturnal_insects; not_day_only_control",source="McMullen 2012 pollination; PMC3489146",night_moth_effective=False,evolutionary_day_to_night_shift_proven=False,endemic_range="not_Galapagos_endemic"),
 "Ipomoea carnea":dict(archipelago="widespread",evidence="both_time_bloom_author_trait; not_effective_day_only; shrub_vs_habeliana_vine",source="Kendall Nicholson 2025 DOI:10.1111/ele.70036",night_moth_effective=False,evolutionary_day_to_night_shift_proven=False,endemic_range="not_archipelago_endemic"),
 "Echium wildpretii":dict(archipelago="Tenerife",evidence="bird_and_insect_visitors; not_strict_day_only",source="doi:10.1016/S1146-609X(02)01167-0",night_moth_effective=False,evolutionary_day_to_night_shift_proven=False,endemic_range="confirmed_Canary_archipelago"),
 "Ceropegia fusca":dict(archipelago="Canary_Islands",evidence="genus_has_fly_pollination; focal_species_direct_effectiveness_not_verified",source="PMCID:PMC2701756",night_moth_effective=False,evolutionary_day_to_night_shift_proven=False,endemic_range="requires_species_level_review"),
 "Roussea simplex":dict(archipelago="Mauritius",evidence="field_experiment_2009_effective_endemic_day_gecko_Phelsuma_cepediana; later_bird_effect_reported; not_an_exclusive_permanent_pollinator_assumption",source="doi:10.1086/593050; Kew POWO Roussea simplex",night_moth_effective=False,gecko_effective=True,evolutionary_day_to_night_shift_proven=False,endemic_range="confirmed_Mauritius_only"),
 "Ruizia blackburniana":dict(archipelago="Mauritius",evidence="effective_gecko_Phelsuma_cepediana_pollination_and_Pandanus_spatial_context; original_study_name_Trochetia_blackburniana_is_WCVP_synonym",source="PMID:17262697; Kew POWO Ruizia blackburniana accepted",night_moth_effective=False,gecko_effective=True,evolutionary_day_to_night_shift_proven=False,endemic_range="confirmed_Mauritius_only"),
}
def require(ok,reason):
    if not ok:raise ValueError("ALAN_ISLAND_NIGHTPAIR: "+str(reason))
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(4<<20),b""):h.update(part)
    return h.hexdigest()
def blob(path):
    p=Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(p)}\0".encode()+p).hexdigest()
def sources(study,traits,pcnl,group_dir):
    from scripts.analysis.alan_diel_exclusion_experimental_source_overlap_20261110 import (
        TRAIT_BLOB,normalize_name,
    )
    require(sha(pcnl)==SITE_SHA,"Original source prior-1992-2024 PCNL file changed")
    assert TRAIT_BLOB==STUDY_SOURCE_BLOB
    require(blob(traits)==STUDY_SOURCE_BLOB,"Independently authored plant traits Git object differs")
    site=pd.read_csv(pcnl,dtype={"island":str})
    require(len(site)==85318 and not site.duplicated(KEY).any(),"Source site/year duplicated or changed")
    valid=site.loc[site.historic_ntl_complete.astype(str).str.lower().eq("true")].copy()
    v=valid.groupby(["island","year"],observed=True).pcnl_1992_to_previous_year_mean_log1p.agg(["size","min","max"])
    keep=set(v.loc[v["size"].ge(2)&(v["max"]-v["min"]).ge(.02)].reset_index().island)
    require(len(keep)==562,"Original 562 island source light contrast support changed")
    matched=valid.loc[valid.island.isin(keep)].copy()
    require(len(matched)==77647 and matched.island.nunique()==562,"Original ALAN macro population changed")
    prior=pd.read_csv(study)
    require({"species","has_day_night_exclusion_study","source_matched_WCVP_accepted"}<=set(prior),"Original PR29 study overlap not provided")
    studied=set(prior.loc[prior.has_day_night_exclusion_study &
                prior.source_matched_WCVP_accepted,"species"])
    require(len(studied)==54,"Kendall original exact-WCVP 54 study plants changed")
    trait=pd.read_csv(traits,dtype=str,keep_default_na=False)
    require(len(trait)==142 and {"plant_species","bloom_period_simple","plant_family","life_form","color"}<=set(trait),"Original author trait evidence changed")
    trait["species"]=normalize_name(trait.plant_species)
    require(trait.species.is_unique,"Multiple unscreened original author trait rows")
    # All original study plants PLUS separately sourced reviewed island probes.
    names=sorted(studied|set(SEEDS))
    chunks=[]
    for g in range(16):
        db=Path(group_dir)/f"allflora-island-group-{g:02d}.sqlite"
        rec=db.with_suffix(".json")
        require(db.is_file() and rec.is_file(),"Original 16 source SQL+receipt unavailable")
        obj=json.loads(rec.read_text())
        require(obj.get("group")==g and obj.get("schema")==
                "ALAN_full93m_genuine_WCVP_allflora_moth_evidence_island_source_group_v1" and
                obj["counts"]["completed_verified_original_tiers"]==14,
                "Original photo-unique signed 14-tier groups changed")
        with sqlite3.connect(db) as con:
            for i in range(0,len(names),70):
                group=names[i:i+70]
                sql="SELECT island,latbin,lonbin,year,species,photo_ids,observation_uuids,moth_positive_photo_ids,butterfly_positive_photo_ids FROM flora WHERE species IN ("+",".join("?" for _ in group)+")"
                a=pd.read_sql_query(sql,con,params=group)
                if len(a):chunks.append(a)
    d=pd.concat(chunks,ignore_index=True)
    d["island"]=d.island.astype(str)
    d=d.groupby(KEY+["species"],observed=True)[["photo_ids","observation_uuids","moth_positive_photo_ids","butterfly_positive_photo_ids"]].sum().reset_index()
    require(d.year.isin(YEARS).all() and d.photo_ids.ge(d.observation_uuids).all(),
            "Original photo/UUID or years malformed")
    return d,matched,valid,studied,trait
def run(study,traits,pcnl,groups,out):
    d,full,valid,studied,traits=sources(study,traits,pcnl,groups)
    keys=KEY
    pcnl_keys=full[keys]
    eligible=d.merge(pcnl_keys,on=keys,how="inner",validate="many_to_one")
    cov=d.merge(valid[keys+["pcnl_1992_to_previous_year_mean_log1p"]],on=keys,how="left",validate="many_to_one")
    cov=cov.loc[cov.pcnl_1992_to_previous_year_mean_log1p.notna()]
    light=cov.groupby(["species","island","year"],observed=True).pcnl_1992_to_previous_year_mean_log1p.agg(["size","min","max"]).reset_index()
    light["spatial_light_contrast"]=light["size"].ge(2)&(light["max"]-light["min"]).ge(.02)
    records=[]
    for name in sorted(studied|set(SEEDS)):
        raw=d.loc[d.species.eq(name)]
        enough=eligible.loc[eligible.species.eq(name)]
        q=light.loc[light.species.eq(name)]
        focus=dict(SEEDS.get(name,{}))
        focus.update(species=name,has_independent_day_night_exclusion_study=name in studied,
            island_original_photo_IDs=int(raw.photo_ids.sum()),
            island_original_independent_observation_UUIDs=int(raw.observation_uuids.sum()),
            island_polygon_units=int(raw.island.nunique()),
            island_0p05degree_site_years=len(raw),
            source_562_PCNL_site_years=len(enough),
            source_562_PCNL_photo_IDs=int(enough.photo_ids.sum()),
            source_562_PCNL_islands=int(enough.island.nunique()),
            island_units_sameyear_light_contrast=int(q.loc[q.spatial_light_contrast,"island"].nunique()),
            sameyear_island_light_contrast_opportunities=int(q.spatial_light_contrast.sum()),
            original_GSHHG_island_IDs=";".join(sorted(raw.island.unique())))
        records.append(focus)
    table=pd.DataFrame(records)
    # Co-recorded independent-experiment species pairs are NOT implicitly night-vs-day.
    experimental=d.loc[d.species.isin(studied)]
    pairrows=[]
    for cell,sub in experimental.groupby(keys,observed=True):
        if len(sub)<2:continue
        for a,b in itertools.combinations(sorted(sub.species.unique()),2):
            pairrows.append((a,b,cell[0],cell[3]))
    pair=pd.DataFrame(pairrows,columns=["species_A","species_B","island","year"])
    groups=pair.groupby(["species_A","species_B"],observed=True).agg(
        co_recorded_site_years=("year","size"),co_recorded_islands=("island","nunique")).reset_index()
    tr=traits.set_index("species")
    for label,axis in [("A","species_A"),("B","species_B")]:
        groups["genus_"+label]=groups[axis].str.split().str[0]
        groups["bloom_"+label]=groups[axis].map(tr.bloom_period_simple).fillna("unknown")
        groups["family_"+label]=groups[axis].map(tr.plant_family).fillna("unknown")
        groups["form_"+label]=groups[axis].map(tr.life_form).fillna("unknown")
    groups["same_genus"]=groups.genus_A.eq(groups.genus_B)
    groups["strict_day_vs_night_bloom"]=(
       (groups.bloom_A.eq("night") & groups.bloom_B.eq("day"))|
       (groups.bloom_B.eq("night") & groups.bloom_A.eq("day")))
    groups["same_genus_night_day"]=groups.same_genus & groups.strict_day_vs_night_bloom
    require(len(groups)==298 and groups.same_genus.sum()==8 and
            groups.same_genus_night_day.sum()==0,
            "Original author trait × exact photo site pairs changed")
    ip=table.loc[table.species.eq("Ipomoea habeliana")].iloc[0]
    require((int(ip.island_original_photo_IDs),int(ip.island_original_independent_observation_UUIDs),
             int(ip.island_polygon_units),int(ip.island_0p05degree_site_years),
             int(ip.source_562_PCNL_site_years),
             int(ip.island_units_sameyear_light_contrast))==(44,22,5,10,3,0),
            "Original island hawkmoth-effective native candidate photo opportunity changed")
    br=table.loc[table.species.eq("Brighamia insignis")].iloc[0]
    require(int(br.island_original_photo_IDs)==0,
            "Hawaiian putative moth-syndrome no-source data unexpectedly changed")
    cc=groups.loc[((groups.species_A=="Cordia lutea") & (groups.species_B=="Ipomoea habeliana")) |
                  ((groups.species_A=="Ipomoea habeliana") & (groups.species_B=="Cordia lutea"))]
    require(len(cc)==1 and int(cc.co_recorded_site_years.iloc[0])==8 and
            int(cc.co_recorded_islands.iloc[0])==3,
            "Galapagos mixed nocturnal control co-recording support not conserved")
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite signed source opportunity audit")
    out.mkdir(parents=True)
    table.to_csv(out/"reviewed_island_night_pollination_seed_photo_support.csv",index=False)
    groups.to_csv(out/"source_exact_site_coobserved_diel_experiment_pairs.csv",index=False)
    report=dict(schema=SCHEMA,source_photo_registry_IDs=93209351,
      original_accepted_experimental_species=54,
      original_experimental_coobserved_species_pairs=len(groups),
      coobserved_same_genus_pairs=int(groups.same_genus.sum()),
      coobserved_same_genus_strict_night_day_bloom_pairs=int(groups.same_genus_night_day.sum()),
      ancestral_day_to_night_shift_verified_for_seed_species=0,
      effective_hawkmoth_island_endemic_covered_I_habeliana_photos=44,
      effective_hawkmoth_island_endemic_covered_I_habeliana_observations=22,
      I_habeliana_562_macro_source_sites=3,
      I_habeliana_sameyear_light_contrast_islands=0,
      Brighamia_insignis_photo_source_IDs=0,
      I_habeliana_plus_Cordia_lutea_coobserved_site_years=8,
      Cordia_lutea_not_a_strict_diurnal_control=True,
      source_group_island_native_status_unverified=True,
      bloom_period_not_effective_pollination_dependency=True,
      directional_island_night_evolution_unverified=True,
      island_vs_mainland_3way_interaction_estimable=False,
      seed_register_not_globally_exhaustive=True)
    (out/"ALAN_actual_island_nightshift_pair_support_gate.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_SOURCE_NIGHT_SHIFT_MATCHED_PAIR_GATE",
          json.dumps(report,sort_keys=True),flush=True)
    return report
def main():
    p=argparse.ArgumentParser()
    for k in ("study","traits","pcnl","groups","out"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    run(a.study,a.traits,a.pcnl,a.groups,a.out)
if __name__=="__main__":main()
