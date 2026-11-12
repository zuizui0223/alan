"""Reconstruct the dated photographs of all 544 shared plants in the frozen cells."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from alan_gecko_primary_fiveyear import primary
from alan_globi_v011_same_fiveyear import ID, YEARS, read_one

ROOT = Path(__file__).resolve().parents[2]


def frame(signed, sites, folder):
    folder.mkdir(exist_ok=True, parents=True)
    cols = ["photo_id","observation_id","inat_taxon_id","observed_on","latitude","longitude","source_subarm"]
    p = read_one(signed/"canonical","canonical*.gz",usecols=cols)
    p = p[p.source_subarm.eq("historical_main") & p.observed_on.str[:4].isin(map(str,YEARS))
          & p.latitude.between(-89.999,89.999) & p.longitude.between(-179.999,179.999)].copy()
    p["key"] = (p.observed_on.str[:4]+"|"+np.floor((p.longitude+180)*20).astype(int).astype(str)
                +"_"+np.floor((p.latitude+90)*20).astype(int).astype(str))
    p = p.merge(sites[["key","region5","land_class","region_year","viirs_median","year"]],
                on="key",validate="m:1").drop(columns=["latitude","longitude","source_subarm"])
    assert p.key.nunique()==10698 and len(p)==sites.photo_attempts.sum()
    threshold = sites.groupby("region_year").viirs_median.median()
    p["light"] = np.where(p.viirs_median.gt(p.region_year.map(threshold)),"bright","dark")
    p["month"] = pd.to_datetime(p.observed_on,format="%Y-%m-%d").dt.month
    p.to_parquet(folder/"universe.parquet")
    tax = read_one(signed/"taxonomy","*resolution.csv.gz",dtype={ID:str})
    names = read_one(signed/"taxonomy","*distribution*.gz",usecols=[ID,"wcvp_accepted_name"],dtype={ID:str}).drop_duplicates()
    shared = pd.read_csv(ROOT/"results/moth_pollen_evidence/shared_544_evidence.csv").plant
    p = p.merge(tax[tax.wcvp_match_state.eq("unique_exact_name_to_accepted_id")][["inat_taxon_id",ID]],
                on="inat_taxon_id",validate="m:1").merge(names,on=ID,validate="m:1")
    p = p[p.wcvp_accepted_name.isin(shared)].rename(columns={"wcvp_accepted_name":"plant"})
    assert set(p.plant)==set(shared) and len(shared)==544
    p.to_parquet(folder/"target.parquet")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("inputs","signed","private_cache"):
        parser.add_argument(key,type=Path)
    args=parser.parse_args()
    sites, _, _, _ = primary(args.inputs,args.signed)
    frame(args.signed,sites,args.private_cache)
