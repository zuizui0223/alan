"""Screen NEON phenology sample locations and the official open teaching subset."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import zipfile

import geopandas as gpd
from mois_global_2020_GSHHG_scope_20261107 import assign_points, prepare_land

import numpy as np
import pandas as pd
import requests
import time

API = "https://data.neonscience.org/api/v0"
ARCHIVE = "https://ndownloader.figshare.com/files/22775042"
ARCHIVE_SHA = "a879e2f9ae6f7fed996fe678f67c5625f858530c84c390ecbca771f6ef701027"
YEARS = {2019, 2021, 2023, 2024, 2025}
OUT = Path(__file__).resolve().parents[2] / "results/phenology_monitoring/neon"


def json_data(url, cache):
    name = cache / (sha256(url.encode()).hexdigest() + ".json")
    if not name.exists():
        for attempt in range(4):
            resp = requests.get(url, timeout=60)
            if resp.status_code != 429 or attempt == 3:
                break
            time.sleep(float(resp.headers.get("Retry-After", 65)))
        resp.raise_for_status()
        name.write_bytes(resp.content)
    return json.loads(name.read_bytes())["data"]


def key(lat, lon):
    return str(int(np.floor((lon + 180) * 20))) + "_" + str(int(np.floor((lat + 90) * 20)))


def run(cache, frame, v011):
    cache.mkdir(exist_ok=True, parents=True)
    OUT.mkdir(exist_ok=True, parents=True)
    u = pd.read_parquet(frame / "universe.parquet").drop_duplicates("key")
    taxa = set(pd.read_parquet(frame / "target.parquet").plant)
    assert len(u) == 10698 and len(taxa) == 544
    get = lambda url: json_data(url, cache)
    product = get(API + "/products/DP1.10055.001")
    sites = {s["siteCode"]: set(int(m[:4]) for m in s["availableMonths"]) & YEARS for s in product["siteCodes"]}
    def locations(addr):
        return get(API + "/locations/" + addr)
    with ThreadPoolExecutor(6) as pool:
        nodes = list(pool.map(locations, sites))
        plot_ids = sorted({p for n in nodes for p in n["locationChildren"] if p.lower().endswith(".phe")})
        plots = list(pool.map(locations, plot_ids))
        regions = set(u.region5)
        def near_region(n):
            lat, lon = n["locationDecimalLatitude"], n["locationDecimalLongitude"]
            return lat is not None and lon is not None and any(
                f"{int(np.floor((lon + dx + 180)/5))}_{int(np.floor((lat + dy + 90)/5))}" in regions
                for dx in (-.05,.05) for dy in (-.05,.05))
        children = sorted({p for n in plots if near_region(n) for p in n["locationChildren"]})
        corners = list(pool.map(locations, children))
    points = []
    regs = set(u.region_year)
    lookup = u.set_index("key")
    for n in plots + corners:
        lat, lon = n["locationDecimalLatitude"], n["locationDecimalLongitude"]
        if lat is None or lon is None:
            continue
        region = str(int(np.floor((lon + 180) / 5))) + "_" + str(int(np.floor((lat + 90) / 5)))
        for year in sorted(sites[n["siteCode"]]):
            grid = str(year) + "|" + key(lat, lon)
            frame_region = str(year) + "|" + region in regs
            if frame_region:
                match = lookup.loc[grid] if grid in lookup.index else None
                points.append((n["siteCode"], n["locationName"], year, region,
                               match.land_class if match is not None else "",
                               match.light if match is not None else ""))
    points = pd.DataFrame(points, columns=["site", "named_location", "year", "region5", "land_class", "light"])
    coords = pd.DataFrame([(n["locationName"], n["locationDecimalLatitude"], n["locationDecimalLongitude"])
                           for n in plots + corners], columns=["unit_id","latitude","longitude"])
    coast = v011 / "coast"
    receipt = json.loads((coast / "source_receipt.json").read_text())
    for name, digest in receipt["extracted_sha256"].items():
        assert sha256((coast / name).read_bytes()).hexdigest() == digest
    assigned = assign_points(coords, prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp")),
                             gpd.read_file(coast/"GSHHS_h_L2.shp"))
    points = points.merge(assigned[["unit_id","land_class"]].rename(columns={"unit_id":"named_location","land_class":"position_land_class"}),on="named_location",validate="m:1")
    points.to_csv(OUT/"plot_screen.csv",index=False)
    archive = cache / "neon_tutorial.zip"
    if not archive.exists():
        response = requests.get(ARCHIVE,timeout=100)
        response.raise_for_status()
        archive.write_bytes(response.content)
    digest = sha256(archive.read_bytes()).hexdigest()
    if ARCHIVE_SHA and ARCHIVE_SHA != digest:
        raise ValueError("Unexpected teaching archive digest")
    obs, indiv = [], []
    with zipfile.ZipFile(archive) as package:
        for name in package.namelist():
            if ".DP1.10055.001." not in name or not name.endswith(".zip") or not any("."+str(y)+"-" in name for y in YEARS):
                continue
            with zipfile.ZipFile(BytesIO(package.read(name))) as month:
                for file in month.namelist():
                    if "phe_perindividual." in file:
                        indiv.append(pd.read_csv(month.open(file)))
                    if "phe_statusintensity." in file:
                        obs.append(pd.read_csv(month.open(file)))
    obs = pd.concat(obs).drop_duplicates().sort_values("editedDate", na_position="first")
    obs = obs[obs.phenophaseName.eq("Open flowers")].drop_duplicates(["individualID","date","phenophaseName"],keep="last")
    obs = obs[obs.phenophaseStatus.isin(["yes","no"])].copy()
    indiv = pd.concat(indiv).drop_duplicates().sort_values("editedDate", na_position="first").drop_duplicates("individualID",keep="last")
    indiv["plant"] = indiv.scientificName.str.extract(r"^([A-Z][a-z]+ [a-z][a-z-]+)")
    wcvp_zip = v011/"wcvp.zip"
    assert sha256(wcvp_zip.read_bytes()).hexdigest() == "d32ea2b3a85e489b14e83bcc9eae7274532e1d113753f7be290d4b2dfde573fa"
    with zipfile.ZipFile(wcvp_zip) as wcvp:
        w = pd.read_csv(wcvp.open("wcvp_names.csv"),sep="|",dtype=str,keep_default_na=False,
                        usecols=["taxon_name","taxon_authors","taxon_rank","taxon_status","plant_name_id","accepted_plant_name_id"])
    w = w[w.taxon_rank.eq("Species")]
    accepted = w[w.taxon_status.eq("Accepted")].set_index("plant_name_id").taxon_name
    mapping = {}
    for name, full_name in indiv[["plant","scientificName"]].drop_duplicates().itertuples(index=False,name=None):
        group = w[w.taxon_name.eq(name)]
        if not isinstance(full_name,str) or not isinstance(name,str):
            continue
        author = full_name[len(name):].replace(" ","")
        authored = group[group.taxon_authors.str.replace(" ","",regex=False).eq(author)]
        if len(authored):
            group = authored
        ids = group.accepted_plant_name_id.unique()
        if len(ids)==1 and ids[0] in accepted and group.taxon_status.isin(["Accepted","Synonym","Homotypic Synonym","Heterotypic Synonym"]).all():
            mapping[full_name] = accepted[ids[0]]
    indiv["plant"] = indiv.scientificName.map(mapping)
    indiv = indiv[indiv.plant.isin(taxa)]
    joint = obs.merge(indiv[["individualID","plant","sampleLatitude","sampleLongitude","decimalLatitude","decimalLongitude"]],on="individualID")
    joint["year"] = joint.date.str[:4].astype(int)
    # Polygon assignment is conservative if a tagged individual lacks sample coordinates.
    fallback = joint.sampleLatitude.isna() | joint.sampleLongitude.isna()
    joint["plot_fallback"] = fallback
    lat = joint.sampleLatitude.where(~fallback, joint.decimalLatitude)
    lon = joint.sampleLongitude.where(~fallback, joint.decimalLongitude)
    joint["key"] = [str(y) + "|" + key(a,b) for y,a,b in zip(joint.year,lat,lon)]
    regional = joint.merge(points[["named_location","year","region5","position_land_class"]].drop_duplicates(),
                           left_on=["namedLocation","year"],right_on=["named_location","year"],validate="m:1")
    reg_obs = regional.groupby(["plant","year","region5","position_land_class","siteID","phenophaseStatus"]).date.agg(["size","nunique"]).reset_index()
    reg_obs.to_csv(OUT/"teaching_region_status.csv",index=False)
    hits = joint.merge(u[["key","region5","land_class","light"]], on="key", validate="m:1")
    if len(hits):
        hits.groupby(["plant","region5","year","land_class","light","plot_fallback","phenophaseStatus"]).date.agg(["size","nunique"]).reset_index().to_csv(OUT/"teaching_arm_status.csv",index=False)
    else:
        pd.DataFrame(columns=["plant","region5","year","land_class","light","plot_fallback","phenophaseStatus","size","nunique"]).to_csv(OUT/"teaching_arm_status.csv",index=False)
    p = points[points.land_class.ne("")]
    summary = {"product":"DP1.10055.001","source_sites":len(sites),"public_phenology_plots":len(plots),"public_plot_child_locations":len(corners),
            "sites_within_frozen_region_year":sorted(points.site.unique()),"regions_within_frozen_region_year":sorted(points.region5.unique()),
            "region_plot_land_classes":points.position_land_class.value_counts().to_dict(),
            "sites_with_point_in_frozen_cell_year":sorted(p.site.unique()),"arms_of_plot_locations":sorted(set(p.land_class+"_"+p.light)),
            "teaching_archive_sha256":digest,"teaching_years":sorted(joint.year.unique().tolist()),"teaching_sites":sorted(obs.siteID.unique()),
            "teaching_yes_no_records":len(obs),"teaching_shared_species":joint.plant.nunique(),"teaching_shared_species_records":len(joint),
            "teaching_cell_year_matched":len(hits),"teaching_cell_year_matched_exact_sample_coordinate":int((~hits.plot_fallback).sum()),
            "teaching_cell_year_matched_species":hits.plant.nunique(),"teaching_cell_year_matched_regions":hits.region5.nunique()}
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(summary)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("frame", type=Path)
    parser.add_argument("v011", type=Path)
    args=parser.parse_args()
    run(args.cache,args.frame,args.v011)
