"""Trace named adult prey to their recorded flower plants in the fixed five-year cohort."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from alan_gecko_primary_fiveyear import ID, MULTIGUILD, positive_records, primary, prepare_land, assign_points
from alan_gecko_moth_species_edges import ARTHROPODS, SOURCE_SHA, SOURCE_ROWS

SCALES = {"1deg": ["x1", "y1"], "0.05deg": ["gx", "gy"]}
BASE = ["year", "land_class", "land_id"]


def prepare(inputs, signed, predators, private):
    private.mkdir(exist_ok=True, parents=True)
    s, digest, plants, legacy = primary(inputs, signed)
    parts = s.target_id.str.split("|", regex=False, expand=True)
    s["land_id"] = parts[2]
    grid = parts[1].str.split("_", expand=True).astype(int)
    s["gx"], s["gy"] = grid[0], grid[1]
    s["x1"], s["y1"] = s.gx // 20 - 180, s.gy // 20 - 90
    needed = {str(y): g[["x1","y1"]].drop_duplicates().values.tolist() for y,g in s.groupby("year")}
    taxa = {k:v for k,v in MULTIGUILD.items() if v in ARTHROPODS}
    land, provenance = positive_records(predators, inputs / "coast", taxa, needed)
    seen = []
    for scale,dims in SCALES.items():
        cols = BASE + dims
        det = land[cols+["species"]].drop_duplicates().rename(columns={"species":"predator"})
        seen.append(s.merge(det,on=cols).assign(scale=scale))
    pd.concat(seen)[["key",*BASE,"gx","gy","x1","y1","predator","scale"]].to_parquet(private/"predator_sites.parquet",index=False)
    plants[["key", ID, "wcvp_accepted_name"]].to_parquet(private/"plants.parquet",index=False)
    visits = legacy[legacy.guild.eq("moth_family_visit")][["insect_name","plant_name"]].drop_duplicates()
    visits.rename(columns={"insect_name":"moth","plant_name":"wcvp_accepted_name"}).to_parquet(private/"visits.parquet",index=False)
    (private/"primary.json").write_text(json.dumps(dict(target_hash=digest, cohort=s.groupby("land_class").size().to_dict(),
                                                     predator_snapshots=provenance["input_sha256"])), encoding="utf8")


def candidates(private, edges):
    receipt = json.loads((edges/"receipt.json").read_text())
    assert receipt["source_sha256"] == SOURCE_SHA and receipt["scanned"] == SOURCE_ROWS
    assert hashlib.sha256((edges/"species_prey.csv").read_bytes()).hexdigest() == receipt["prey_table_sha256"]
    prey = pd.read_csv(edges/"species_prey.csv", keep_default_na=False)
    assert set(prey.predator) == set(ARTHROPODS)
    adult = prey.loc[prey.source_stage.eq("adult") & prey.moth.ne(""),["predator","moth"]].drop_duplicates()
    visits = pd.read_parquet(private/"visits.parquet")
    links = adult.merge(visits,on="moth")
    plants = pd.read_parquet(private/"plants.parquet")
    sites = pd.read_parquet(private/"predator_sites.parquet")
    present = links.merge(plants,on="wcvp_accepted_name")
    possible = present.merge(sites,on=["key","predator"])
    return adult, links, present, possible


def summarize(matches):
    return {group: int(matches.loc[matches.land_class.eq(group),"key"].nunique()) for group in ("island","mainland")}


def finish(inputs, private, edges, output):
    adult, links, present, possible = candidates(private,edges)
    receipts = json.loads((private/"moth_queries.json").read_text())
    needed = possible.query("scale=='1deg'")[["moth","year","x1","y1"]].drop_duplicates()
    checked = set()
    records = {}
    for q in receipts:
        with gzip.open(private/"moth_queries"/q["file"],"rb") as f:
            raw = f.read()
        assert hashlib.sha256(raw).hexdigest() == q["sha256"]
        obj = json.loads(raw)
        assert len(obj["rows"]) == q["count"]
        for x,y in q["cells"]:
            checked.add((q["moth"],q["year"],x,y))
        for r in obj["rows"]:
            assert r["year"] == q["year"] and r["basisOfRecord"]=="HUMAN_OBSERVATION"
            if r["species"] != q["moth"] or r["speciesKey"] != q["taxon"]: continue
            if r["occurrenceStatus"] not in ("PRESENT","DETECTED") or not 0 < (r["coordinateUncertaintyInMeters"] or 0) <= 1000: continue
            if str(r["lifeStage"]).lower() not in ("adult","imago"): continue
            records[r["key"]] = r
    assert checked == set(needed.itertuples(index=False,name=None))
    datasets = sorted({r["datasetKey"] for r in records.values() if r["datasetKey"]})
    d = pd.DataFrame(records.values()).rename(columns={"key":"unit_id","species":"moth","decimalLongitude":"longitude","decimalLatitude":"latitude"})
    if len(d):
        coast = inputs/"coast"
        ass = assign_points(d,prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp")),gpd.read_file(coast/"GSHHS_h_L2.shp"))
        d = d.merge(ass,on="unit_id",validate="one_to_one")
        d = d[d.land_class.isin(["island","mainland"])].copy()
        d["gx"],d["gy"] = np.floor((d.longitude+180)*20).astype(int), np.floor((d.latitude+90)*20).astype(int)
        d["x1"],d["y1"] = np.floor(d.longitude).astype(int),np.floor(d.latitude).astype(int)
    else:
        d=pd.DataFrame(columns=["moth",*BASE,"gx","gy","x1","y1"])
    matched=[]
    for scale,dims in SCALES.items():
        cols=["moth",*BASE,*dims]
        matched.append(possible[possible.scale.eq(scale)].merge(d[cols].drop_duplicates(),on=cols,validate="many_to_one"))
    m=pd.concat(matched,ignore_index=True)
    identity=["predator","moth","key",ID]
    fine=set(m.loc[m.scale.eq("0.05deg"),identity].itertuples(index=False,name=None))
    coarse=set(m.loc[m.scale.eq("1deg"),identity].itertuples(index=False,name=None))
    assert fine <= coarse
    receipt=json.loads((edges/"receipt.json").read_text())
    assert hashlib.sha256((edges/"predator_flowers.csv").read_bytes()).hexdigest()==receipt["predator_flower_sha256"]
    flowers=pd.read_csv(edges/"predator_flowers.csv")
    stages = pd.read_csv(edges/"species_prey.csv", keep_default_na=False)
    totals=[]
    pairs=[]
    for predator in ARTHROPODS:
        group=m[m.predator.eq(predator)]
        f=flowers[flowers.predator.eq(predator)]
        direct=set(f.loc[f.role.eq("reported_pollinates"),"plant"])
        entry=dict(predator=predator, named_adult_prey=adult.query("predator==@predator").moth.nunique(),
                   prey_with_visit_link=links.query("predator==@predator").moth.nunique(),
                   plant_names=links.query("predator==@predator").wcvp_accepted_name.nunique(),
                   plant_names_on_cohort=present.query("predator==@predator").wcvp_accepted_name.nunique(),
                   predator_pollinates_plant_names=sorted(direct),
                   same_plant_with_pollinates_label=sorted(direct & set(links.query("predator==@predator").wcvp_accepted_name)))
        counted = stages[stages.predator.eq(predator)].groupby("source_stage").relationship_rows.sum()
        entry["prey_rows_by_stage"] = counted.to_dict()
        entry["predator_flower_visit_plant_names"] = sorted(set(f.loc[f.role.eq("flower_visit"),"plant"]))
        for scale in SCALES:
            g=group[group.scale.eq(scale)]
            entry[scale]=dict(site_years=summarize(g), moths=g.moth.nunique(), plants=g[ID].nunique(),
                              same_plant_reported_pollinates=summarize(g[g.wcvp_accepted_name.isin(direct)]))
        totals.append(entry)
        for name in sorted(set(adult.query("predator==@predator").moth)):
            v=links[(links.predator==predator)&(links.moth==name)]
            row=dict(predator=predator,moth=name,plant_names=sorted(set(v.wcvp_accepted_name)))
            for scale in SCALES:
                subset = group[(group.moth==name)&(group.scale==scale)]
                row[scale]=summarize(subset)
                if scale == "0.05deg":
                    row["plants_at_matched_sites"] = {side: sorted(set(subset.loc[subset.land_class.eq(side),"wcvp_accepted_name"])) for side in ("island","mainland")}
            pairs.append(row)
    covered={scale:summarize(m[m.scale.eq(scale)]) for scale in SCALES}
    out=dict(primary=json.loads((private/"primary.json").read_text()), source=receipt, totals=totals, union_of_three=covered,
             pairs=pairs, eligible_gbif_dataset_keys=datasets, gbif= [{k:v for k,v in q.items() if k not in ("cells","file")} for q in receipts])
    output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(totals,covered)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest="command",required=True)
    for command,args in (("prepare",["inputs","signed","predators","private"]),("finish",["inputs","private","edges","output"])):
        sub=subs.add_parser(command)
        for arg in args: sub.add_argument(arg,type=Path)
    args=vars(parser.parse_args())
    {"prepare":prepare,"finish":finish}[args.pop("command")](**args)
