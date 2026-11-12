"""Assemble explicit status visits and audit repeated individuals in frozen 0.05-degree cells."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from alan_globi_v011_same_fiveyear import YEARS

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/phenology_repeated_npn"


def collapse(folder):
    expected = [folder / f"open_{year}_{month:02}.parquet" for year in YEARS for month in range(1,13)]
    source = pd.concat([pd.read_parquet(p) for p in expected], ignore_index=True)
    assert source.observation_id.is_unique
    names = pd.read_csv(folder / "target_names.csv")
    source = source.merge(names[["species_id","accepted_name"]],on="species_id",validate="m:1")
    assert source.observation_id.notna().all()
    known = source.loc[source.phenophase_status.isin([0,1])].rename(columns={
        "site_id":"station", "individual_id":"individual","accepted_name":"plant",
        "observation_date":"date","phenophase_status":"status"})
    known = known[known.individual.gt(0) & known.station.gt(0) & known.latitude.between(-90,90) & known.longitude.between(-180,180)]
    group = ["station","individual","plant","date"]
    counts = known.groupby(group).agg(statuses=("status","nunique"), latitudes=("latitude","nunique"),longitudes=("longitude","nunique"))
    agreed = counts[counts.eq(1).all(axis=1)].reset_index()[group]
    visits = known.merge(agreed,on=group,validate="m:1").drop_duplicates(group)
    visits["year"] = visits.date.str[:4].astype(int)
    # Intensity is the proportion of fresh floral structures that are open.
    # Separate possible strength of flowering from the yes/no state.
    cat = known.loc[known.status.eq(1) & known.intensity_value.isin(
        ["Less than 5%", "5-24%", "25-49%", "50-74%", "75-94%", "95% or more"]), group+["intensity_value"]].drop_duplicates()
    cat = cat[~cat.duplicated(group, keep=False)]
    visits = visits[["year","latitude","longitude","plant","station","individual","date","status"]].merge(cat,on=group,how="left",validate="1:1")
    visits.to_parquet(folder / "individual_visits.parquet",index=False)
    names.to_csv(OUT / "species_crosswalk.csv",index=False)
    return visits, dict(catalog_matched_names=int(names.accepted_name.nunique()), records=len(source),
        unknown_status_records=int((~source.phenophase_status.isin([0,1])).sum()),
        conflicting_same_day_or_coordinate_groups=int((~counts.eq(1).all(axis=1)).sum()),
        informative_visits=len(visits), nationwide_plants=int(visits.plant.nunique()),
        files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in expected})


def analyze(visits, universe, summary):
    frame = pd.read_parquet(universe,columns=["key","region5","land_class","light","viirs_median"]).drop_duplicates()
    assert len(frame)==10698 and frame.key.is_unique
    v = visits.assign(key=visits.year.astype(str)+"|"+np.floor((visits.longitude+180)*20).astype(int).astype(str)
                      +"_"+np.floor((visits.latitude+90)*20).astype(int).astype(str)).merge(frame,on="key",validate="m:1")
    ids = ["region5","year","land_class","light","plant","station","individual"]
    v["when_no"] = v.date.where(v.status.eq(0),None)
    v["when_yes"] = v.date.where(v.status.eq(1),None)
    # Dates are ISO local calendar strings. Missing status-specific dates must remain missing.
    for col in ["when_no","when_yes"]:
        v[col] = pd.to_datetime(v[col])
    units = v.groupby(ids).agg(n_days=("date","nunique"),n_yes=("when_yes","count"),n_no=("when_no","count"),
        first_yes=("when_yes","min"),last_yes=("when_yes","max"),first_no=("when_no","min"),last_no=("when_no","max")).reset_index()
    units["repeated_dates"] = units.n_days.ge(2)
    units["both_statuses"] = units.n_yes.gt(0) & units.n_no.gt(0)
    units["bracketed"] = (units.first_no<units.first_yes) & (units.last_no>units.last_yes)
    arms = units.groupby(ids[:5]).agg(individuals=("individual","size"),repeated_dates=("repeated_dates","sum"),
        both_statuses=("both_statuses","sum"),bracketed=("bracketed","sum")).reset_index()
    arms.to_csv(OUT / "exact_cell_arms.csv",index=False)
    counts = {}
    for label in ["individuals","repeated_dates","both_statuses","bracketed"]:
        subset = arms[arms[label]>0]
        all_four = subset.groupby(["plant","region5","year"]).size()
        complete = all_four[all_four.eq(4)].reset_index(name="arms")
        counts[label] = dict(species=int(subset.plant.nunique()),regions=int(subset.region5.nunique()),
            complete_species=int(complete.plant.nunique()), complete_regions=int(complete.region5.nunique()),
            complete_strata=len(complete))
    summary.update(exact_cell_informative_visits=len(v), exact_cell_individual_site_years=len(units),exact_cells=v.key.nunique(),coverage=counts)
    (OUT / "summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({k:t for k,t in summary.items() if k!="files_sha256"},indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("private_cache",type=Path)
    parser.add_argument("universe",type=Path)
    args=parser.parse_args()
    OUT.mkdir(exist_ok=True,parents=True)
    visits,summary=collapse(args.private_cache)
    analyze(visits,args.universe,summary)
