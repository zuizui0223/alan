#!/usr/bin/env python3
"""ALAN independently new observation UUIDs from the immutable full93m tier06.

Source is the *original* 597,070-photo observation-hash tier, selected before
any colour ROI. Its 16 asset SHA256 values and global registry are pinned.
Do not download all 93m images, infer colours, fabricate flowering/native
status, or add multiple photographs of the same biological observation.
This is an ALAN-specific conditional acquisition queue, not a random
global-flower sample. Existing five-year models remain unchanged.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,re,tarfile
from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import prepare_land,assign_points
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import cell_id
from scripts.analysis.mois_global2025_original_source_20261107 import COAST_SHA

REGISTRY_SHA="4ed5139e510ba8442d4cdb03f2b769d94fac942832b1d544d7a7925bc431255b"
TAR_SHA=(
"79553cd08d82387653815e0df568ed21bbc7514ae4b3de2466e8e6d67bae3211",
"e6cb5e1d5436370c6b883b610246f33a7a95c81f5940a1847b208f54107a2487",
"db34f2f7de0f035e923d52639c139edffedd3f1b70749de0d64355e3917b02f7",
"a2d69ff864da0244f18ceae6ca4cc1e7f98c5121f26b149f1e7ca9fb213e323d",
"a6bdfc9def9decf5e730405c6a2a19cd0d04a1e519324f84bc0e8c439dc83fc2",
"c580fbc0495ea9a0b2c91efac2cf29092aab7f21fde1ecfe5f88b816df9d3f52",
"8787b61ea54146b888c1e3a572bec574c9ed028794718184cab19269255131df",
"02ea26db2f2c1386bf5825b3d493b1d1365c2859eaf1056a9d6f0ab8917239e",
"450c7602a7ecba9c0caf1e9fc8e8c26fc7a7662cce0f423fd20c201ad729ea8b",
"029eab365fa2479abcee8f5e975e772cf51f6eef2d77aa1a839b1f57479c2658",
"8b88141a61eccf2534e75ba98c7fdfa91084baa90d9cbad0de31da062dd4e47e",
"6027c1b3edeb54375617d39cc07f72a67cbd0455518ee1307fe2f2d5e95bdfd3",
"5b99039f23fb7eb65d8824f4fda71f5262bdf7026bcac7bc46f6399df9161cc8",
"1876bb02722fce802f8ef5b4bba62afe91d73104b33288456f52dc106edb3ea26",
"8abd206a48d85318414bf15242635ee35bd32ba715fb694b6207eb4b20404215",
"38a1dabe94bc7cfe9be5c272424ce6067c534b028c15d1a618fcaaf1")
COLUMNS=("photo_id","observation_uuid","inat_taxon_id","observer_id",
"photographer_id","observed_on","latitude","longitude","positional_accuracy",
"photo_license","image_extension","frozen_loader_supported","valid_observed_date",
"valid_geopoint","observation_inclusion_score","tier")
YEARS=(2019,2021,2023,2024,2025)
PHOTO_SALT="ALAN_93m_original_new_obs_photo_blind_v1"
ARM="alan_93m_tier06_new_island_observation_ROI_source_v1"
QUEUE_SIZE=256

def choose_per_observation(rows, known, priority_salt=PHOTO_SALT):
    """Deterministic, colour-blind one-photo selection for each NEW UUID."""
    selected={}
    previous=set(map(str,known))
    for row in rows:
        obs=str(row["observation_uuid"])
        if obs in previous:continue
        pid=str(row["photo_id"])
        score=hashlib.sha256(f"{priority_salt}|{obs}|{pid}".encode()).hexdigest()
        if obs not in selected or score<selected[obs][0]:
            selected[obs]=(score,dict(row))
    return [x[1] for x in selected.values()]

def read_tier(directory,receipt_file,known):
    panel.require(panel.sha(receipt_file)==REGISTRY_SHA,
                  "93m source full-tier registry SHA changed")
    receipt=json.loads(Path(receipt_file).read_text())
    panel.require(receipt.get("schema")==
                  "mois_complete_full93m_original_photo_target_registry_v1"
                  and receipt.get("exact_original_photo_IDs_in_all_disjoint_tiers")==93209351,
                  "Not the published 93m source-only registry")
    all_rows=[];ids=set();counts=[]
    for group,expected_sha in enumerate(TAR_SHA):
        path=Path(directory)/f"tier-06-group-{group:02d}.tar"
        panel.require(path.is_file() and panel.sha(path)==expected_sha,
                      f"Frozen tier06 group {group} bytes modified")
        count=0;parts=set()
        with tarfile.open(path,"r:") as tar:
            for m in tar.getmembers():
                if not m.isfile() or not m.name.endswith(".csv.gz"):continue
                match=re.search(r"part-(\d{3})\.csv\.gz$",m.name)
                panel.require(match is not None,"Non-source photo CSV in tier archive")
                part=int(match.group(1))
                panel.require(part//8==group and part not in parts,
                              "Original taxon-modulo source part wrong")
                parts.add(part)
                with gzip.GzipFile(fileobj=tar.extractfile(m)) as gs:
                    reader=csv.DictReader(io.TextIOWrapper(gs,encoding="utf-8",newline=""))
                    panel.require(tuple(reader.fieldnames or ())==COLUMNS,
                                  "Signed tier source columns changed")
                    for row in reader:
                        pid=row["photo_id"]
                        panel.require(pid.isdecimal() and pid not in ids and row["tier"]=="6",
                                      "Repeated photo or wrong observation-hash tier")
                        ids.add(pid);count+=1
                        all_rows.append(row)
        panel.require(len(parts)==8,"Incomplete source photo subpart group")
        counts.append(count)
    panel.require(len(ids)==len(all_rows)==597070,"Tier06 all sixteen original photo groups incomplete")
    selected=choose_per_observation(all_rows,known)
    return selected,{"upstream_photo_candidates":93209351,
        "tier06_genuinely_indexed_photo_ids":len(all_rows),
        "all_source_tier06_group_counts":counts,
        "novel_original_observation_UUIDs":len(selected),
        "previous_original_canonical_UUIDs_excluded":len(known),
        "one_colour_blind_photo_per_novel_observation":True,
        "source_colour_pixels_or_ROIs":0,
        "source_observation_sampling_probability":0.0064}

def source_gshhg(rows,coast):
    d=pd.DataFrame(rows)
    d["year"]=pd.to_numeric(d.observed_on.astype(str).str.slice(0,4),
                             errors="coerce")
    d["latitude"]=pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"]=pd.to_numeric(d.longitude,errors="coerce")
    valid=(d.year.isin(YEARS)&d.latitude.between(-90,90)&
           d.longitude.between(-180,180))
    basic={"novel_observations_with_selected_year_geocoordinates":int(valid.sum()),
           "all_new_original_observations":len(d)}
    d=d.loc[valid].copy().rename(columns={"photo_id":"unit_id"})
    source=Path(coast)
    rec=json.loads((source/"source_receipt.json").read_text())
    panel.require(rec.get("original_archive_sha256")==COAST_SHA,
                  "GSHHG coastline source archive drift")
    for rel,digest in rec["extracted_sha256"].items():
        panel.require(panel.sha(source/rel)==digest,
                      "GSHHG underlying original polygon SHA drift")
    land=prepare_land(gpd.read_file(source/"GSHHS_h_L1.shp"))
    lakes=gpd.read_file(source/"GSHHS_h_L2.shp")
    pts=assign_points(d,land,lakes)
    d=d.merge(pts[["unit_id","land_id","land_class"]],
              on="unit_id",validate="one_to_one")
    d=d.loc[d.land_class.isin(("island","mainland"))&
              d.land_id.notna()].copy()
    d["land_id"]=d.land_id.astype(str)
    d["year"]=d.year.astype(int)
    d["cell_id"]=cell_id(d.longitude,d.latitude)
    d["target_id"]=d.year.astype(str)+"|"+d.cell_id+"|"+d.land_id
    panel.require(d.unit_id.is_unique and d.observation_uuid.is_unique,
                  "New source repeated photographs or observation UUIDs")
    return d,basic

def environmental_support(source,canonical,wcvp,paths):
    # Reuse existing SHA-gated actual 2019/21/23/24/25 VIIRS/WorldPop numeric
    # site targets. This does NOT assert satellite observations at NEW sites.
    lookup=wcvp.loc[wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
          ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    lookup.inat_taxon_id=lookup.inat_taxon_id.astype(str)
    panel.require(lookup.inat_taxon_id.is_unique,"WCVP accepted ID duplicates")
    years={}
    for year in (2019,2021,2023):
        years[year]=panel.year_historic(year,paths,canonical,wcvp,return_all=True)
    years[2024]=panel.year2024(paths,canonical,wcvp,return_all=True)
    years[2025]=panel.year2025(paths,canonical,wcvp,return_all=True)
    env=[];strict=[]
    for year,(full,existing) in years.items():
        q=full.copy()
        status=q["population_status"] if "population_status" in q else q["population_source_status"]
        density=q["population_density_per_km2"] if "population_density_per_km2" in q else q["population_density"]
        q["environment_valid"]=(q.viirs_valid_months.ge(3)&q.viirs_median.notna()&
            status.eq("OK")&density.notna()&
            q[["bio5","bio12","elevation"]].notna().all(axis=1))
        env.append(q[["target_id","environment_valid"]].copy())
        used=existing[["target_id","land_id","land_class"]].copy()
        used["year"]=year
        strict.append(used)
    e=pd.concat(env,ignore_index=True)
    panel.require(e.target_id.is_unique,"Original env photo cells duplicated")
    q=source.merge(lookup,on="inat_taxon_id",how="left",validate="many_to_one")
    q=q.merge(e,on="target_id",how="left",validate="many_to_one")
    q["has_original_year_ALAN_environment"]=q.environment_valid.fillna(False).astype(bool)
    q["accepted_taxon"]=q.wcvp_accepted_plant_name_id.notna()
    q["supported_image_extension"]=q.image_extension.str.lower().isin(("jpg","jpeg","png"))
    q["source_loader_supported"]=q.frozen_loader_supported.astype(str).str.lower().eq("true")
    r=pd.concat(strict,ignore_index=True)
    original_island_cells=r.loc[r.land_class.eq("island")].groupby("land_id").size()
    q["existing_strict_island_cell_years"]=q.land_id.map(original_island_cells).fillna(0).astype(int)
    return q

def choose_early_island_queue(q,quota=QUEUE_SIZE,max_each=8):
    possible=q.loc[q.land_class.eq("island")&q.accepted_taxon&
        q.has_original_year_ALAN_environment&
        q.supported_image_extension&q.source_loader_supported].copy()
    possible["rank"]=possible.observation_uuid.map(
        lambda s:hashlib.sha256(("ALAN_93m_island_priority_balanced_v1|"+str(s)).encode()).hexdigest())
    possible=possible.sort_values(["existing_strict_island_cell_years","rank"],
                                  kind="stable")
    possible["within_island_rank"]=possible.groupby("land_id").cumcount()
    possible=possible.loc[possible.within_island_rank.lt(max_each)]
    chosen=possible.sort_values(["within_island_rank",
        "existing_strict_island_cell_years","rank"],kind="stable").head(quota).copy()
    panel.require(chosen.unit_id.is_unique and chosen.observation_uuid.is_unique,
                  "New island original acquisition photo/UUID duplicates")
    return chosen

def counts(q):
    return {"original_observations":len(q),
            "distinct_0p05degree_photo_cell_years":int(q.target_id.nunique()),
            "distinct_island_landmasses":int(q.loc[q.land_class.eq("island"),"land_id"].nunique())}

def run(paths,coast,source,registry,out):
    out=Path(out);panel.require(not out.exists(),"Refuse overwrite of fixed pre-colour island queue")
    canonical,wcvp=panel.load_frozen_photos(paths)
    rows,tier_audit=read_tier(source,registry,set(canonical.observation_id.astype(str)))
    geo,geo_audit=source_gshhg(rows,coast)
    ready=environmental_support(geo,canonical,wcvp,paths)
    accepted=ready.accepted_taxon
    island=ready.land_class.eq("island")
    valid=ready.has_original_year_ALAN_environment
    loadable=ready.supported_image_extension&ready.source_loader_supported
    chosen=choose_early_island_queue(ready)
    panel.require(len(chosen)>=64,"Insufficient independent original island observations with actual VIIRS/WorldPop")
    out.mkdir(parents=True)
    queue=chosen.copy()
    queue["photo_id"]=queue.unit_id.astype(str)
    queue["observation_id"]=queue.observation_uuid.astype(str)
    queue["measurement_id"]=queue.photo_id.map(lambda s:"alan93m_"+s)
    queue["source_arm"]="alan93m_new_island_tier06_conditional_priority"
    queue["image_filename"]=queue.measurement_id+".jpg"
    queue["photo_url_large"]=queue.apply(lambda r:
         f"https://inaturalist-open-data.s3.amazonaws.com/photos/{r.photo_id}/large.{r.image_extension.lower()}",axis=1)
    queue["source_snapshot_date"]="2026-09-27"
    queue["selection_probability"]="" # purposive island/coverage queue, not a known π
    queue["shard"]=np.arange(len(queue))//128
    # Original source metadata-only source queue, never flower colours.
    columns=["measurement_id","photo_id","inat_taxon_id","observation_id",
       "observer_id","observed_on","latitude","longitude","positional_accuracy",
       "photo_license","photo_url_large","image_filename","source_arm",
       "selection_probability","source_snapshot_date","shard",
       "land_id","target_id","year","existing_strict_island_cell_years"]
    queue[columns].to_csv(out/"real_new_independent_island_photo_queue.csv.gz",
       compression="gzip",index=False)
    groups=ready.groupby(["year","land_class"],observed=True).agg(
       source_observations=("unit_id","size"),
       accepted_WCVP_taxon_observations=("accepted_taxon","sum"),
       existing_ALAN_environment=("has_original_year_ALAN_environment","sum")
       ).reset_index()
    groups.to_csv(out/"new_original_observation_support_by_year_and_land.csv",index=False)
    summary={
       "schema":"ALAN_93m_original_tier06_256_new_observation_island_source_v1",
       "signed_source_tier06":tier_audit,
       "original_geographic_observation_source":geo_audit,
       "all_novel_source_observations_with_known_land":counts(ready),
       "new_independent_island_observations":counts(ready.loc[island]),
       "new_independent_island_observations_with_accepted_taxon":counts(ready.loc[island&accepted]),
       "new_independent_island_accepted_with_original_year_VIIRS_WorldPop":counts(ready.loc[island&accepted&valid]),
       "new_independent_island_accepted_loadable_with_existing_original_year_environment":counts(ready.loc[island&accepted&valid&loadable]),
       "first_256_blind_new_original_photo_queue":counts(chosen),
       "first_256_source_queue_original_photo_SHA256":panel.sha(out/"real_new_independent_island_photo_queue.csv.gz"),
       "source_photo_per_observation_selected_before_colour":True,
       "source_image_pixels_downloaded":0,"original_photo_ROI_inferences":0,
       "new_species_presence_without_real_photo":0,
       "scope":"separate prioritized acquisition cohort, NEVER append to canonical original historical_main estimand without acquisition-design adjustment",
       "underlying_source_SHA256":{**{k:panel.sha(v) for k,v in paths.items()},
             "coast_original_SHA256":COAST_SHA,
             "93m_original_tier_registry_SHA256":REGISTRY_SHA}
    }
    (out/"ALAN_new_source_island_observation_manifest_receipt.json").write_text(
       json.dumps(summary,sort_keys=True,indent=2)+"\n")
    print("ALAN_GENUINELY_NEW_93M_INDEXED_ISLAND_OBSERVATIONS",
         json.dumps(summary,sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    keys=("canonical","wcvp",
      "geo2019","viirs2019","pop2019",
      "geo2021","viirs2021","pop2021",
      "geo2023","viirs2023","pop2023",
      "geo2024","cn2024","uc2024","uc2024receipt",
      "geo2025","geo2025receipt","viirs2025")
    for key in keys:p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    for arg in ("coast","tier-source","tier-registry","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in keys},a.coast,a.tier_source,a.tier_registry,a.out)
