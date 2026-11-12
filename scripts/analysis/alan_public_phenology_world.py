"""Screen Swedish citizen-science plant calendar and the archived Australian ClimateWatch release."""
import argparse
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "results/public_phenology_world"
YEARS = {2019, 2021, 2023, 2024, 2025}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wcvp_species(zip_path, names):
    columns = ["taxon_name", "taxon_rank", "taxon_status", "plant_name_id", "accepted_plant_name_id"]
    matches = []
    with zipfile.ZipFile(zip_path) as archive, archive.open("wcvp_names.csv") as handle:
        for chunk in pd.read_csv(handle, sep="|", usecols=columns, dtype=str, chunksize=200000):
            matches.append(chunk[chunk.taxon_name.isin(names) & chunk.taxon_rank.eq("Species")])
    taxa = pd.concat(matches)
    taxa = taxa[taxa.taxon_status.isin(["Accepted", "Synonym", "Homotypic Synonym"])]
    taxa["accepted"] = taxa.accepted_plant_name_id.where(taxa.taxon_status.ne("Accepted"), taxa.plant_name_id)
    pairs = taxa[["taxon_name", "accepted"]].dropna().drop_duplicates()
    return pairs[~pairs.taxon_name.duplicated(False)]


def run(folder, wcvp):
    OUT.mkdir(exist_ok=True)
    cohort = pd.read_parquet(folder/"target.parquet")[["plant", "wcvp_accepted_plant_name_id"]].drop_duplicates()
    cells = pd.read_parquet(folder/"universe.parquet")[["key", "region5", "year", "land_class", "light"]].drop_duplicates()
    sw_file = folder/"plants_calendar_2008-2024.csv"
    sw = pd.read_csv(sw_file)
    sw = sw[sw.obs_year.isin(YEARS) & sw.phase_id.isin([51,52,53,54,56])].copy()
    sw["status"] = sw.phase_id.isin([52,53,54]).astype(int)
    sw["date"] = pd.to_datetime(sw.obs_year.astype(str) + sw.obs_doy.astype(str).str.zfill(3), format="%Y%j")
    sw["year"] = sw.obs_year
    sw["region5"] = (np.floor((sw.longitude+180)/5).astype(int).astype(str)+"_"+
                     np.floor((sw.latitude+90)/5).astype(int).astype(str))
    sw["key"] = (sw.year.astype(str)+"|"+np.floor((sw.longitude+180)*20).astype(int).astype(str)+"_"+
                 np.floor((sw.latitude+90)*20).astype(int).astype(str))
    sw = sw.merge(cells[["region5","year"]].drop_duplicates(),on=["region5","year"],validate="m:1")
    raw_count = len(sw)
    names = wcvp_species(wcvp, sw.scientific_name.unique())
    cohort["accepted"] = cohort.wcvp_accepted_plant_name_id.astype(str)
    sw = sw.merge(names,left_on="scientific_name",right_on="taxon_name",validate="m:1").merge(cohort,on="accepted",validate="m:1")
    sw["station"] = sw.latitude.map(str)+"|"+sw.longitude.map(str)
    sw["individual"] = None
    sw["sampling_unit"] = "reported_location_no_individual_id"
    standard_cols = ["year","latitude","longitude","plant","station","individual","date","status","sampling_unit","phase_id","obs_id","quality"]
    sw[standard_cols].to_parquet(folder/"sweden_all_visits.parquet",index=False)
    sw.loc[sw.quality.eq(1),standard_cols].to_parquet(folder/"sweden_standard_visits.parquet",index=False)
    sw.to_parquet(folder/"swedish_phenology_region_year.parquet",index=False)
    frozen = sw.merge(cells,on=["key","region5","year"],validate="m:1")
    def counts(x):
        return {"rows":len(x), "species":x.plant.nunique(), "regions":x.region5.nunique(),
                "region_years":len(x[["region5","year"]].drop_duplicates()), "negative":int(x.status.eq(0).sum()),
                "positive":int(x.status.eq(1).sum())}
    location = ["region5","year","plant","latitude","longitude"]
    repeat = sw[sw.quality.eq(1)].groupby(location).agg(n_days=("date","nunique"), n_states=("status","nunique"))
    both = repeat[(repeat.n_days>=2)&(repeat.n_states==2)].reset_index()
    positive_and_negative = sw[sw.quality.eq(1)].merge(both,on=location,validate="m:1")
    positive_and_negative.to_parquet(folder/"swedish_repeated_region_year.parquet",index=False)
    arm = ["region5","year","plant","land_class","light"]
    frozen.to_parquet(folder/"swedish_phenology_cell_year.parquet",index=False)
    ag = frozen.groupby(arm).agg(records=("date","size"),dates=("date","nunique"),states=("status","nunique")).reset_index()
    ag.to_csv(OUT/"swedish_candidate_arms.csv",index=False)
    archive = folder/"climatewatch.zip"
    with zipfile.ZipFile(archive) as z:
        fields = ET.fromstring(z.read("meta.xml")).find("{*}core").findall("{*}field")
        idx = int(next(f.get("index") for f in fields if f.get("term").endswith("/year")))
        australian = pd.read_csv(z.open("occurrence.tsv"), sep="\t", header=None, usecols=[idx]).iloc[:,0]
    meta = dict(sweden_sha256=sha(sw_file),wcvp_sha256=sha(wcvp),swedish_five_year_flower_rows_in_frozen_region_years=raw_count,
                regional_shared_species=counts(sw), trained_regional_shared_species=counts(sw[sw.quality.eq(1)]),
                regional_repeated_location_both_states_trained=counts(positive_and_negative),same_frozen_cells=counts(frozen),
                australian_archive_sha256=sha(archive),australian_archive_rows=len(australian),
                australian_max_year=int(australian.max()),australian_missing_years=int(australian.isna().sum()),
                australian_focal_year_rows=int(australian.isin(YEARS).sum()))
    (OUT/"data_screen.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(meta)
if __name__=="__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder",type=Path)
    ap.add_argument("wcvp",type=Path)
    args = ap.parse_args()
    run(args.folder,args.wcvp)
