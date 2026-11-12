"""Identify regional open-flower pairs and summarize the extended station time series."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import theilslopes
from alan_long_phenology import OUT


def run(folder):
    plants=pd.read_csv(OUT.parent/"species_interpretation/island_sign_reproductive_evidence.csv").query("island_pp<0")
    visits=pd.read_parquet(folder/"npn/individual_visits.parquet")
    exposure=pd.read_parquet(folder/"npn_exposure.parquet")
    records=visits[visits.status.eq(1)].merge(plants[["plant"]]).merge(exposure,on=["latitude","longitude","year"],validate="m:1")
    by=records.groupby(["region5","plant","year"]).land_class.nunique()
    pairs=by[by.eq(2)].reset_index().groupby(["region5","plant"]).year.agg(["count","min","max"])
    pairs.to_csv(OUT/"frozen_flowering_pair_screen.csv")
    print("regional screen",len(pairs),pairs.index.get_level_values("plant").nunique(),pairs.index.get_level_values("region5").nunique())
    panel=pd.read_csv(OUT/"station_year.csv")
    rows=[]
    for (plant,side),g in panel.groupby(["plant","land"]):
        for field in ["onset_midpoint","ndvi_on_open_dates","viirs","flower_minus_green_days"]:
            valid=g.dropna(subset=[field]).sort_values("year")
            if not len(valid): continue
            result=theilslopes(valid[field],valid.year) if len(valid)>=3 else None
            rows.append(dict(plant=plant,land=side,metric=field,n_years=len(valid),first_year=valid.year.min(),
                last_year=valid.year.max(),median=valid[field].median(),slope_per_year=result.slope if result else None,
                low=result.low_slope if result else None,high=result.high_slope if result else None))
    table=pd.DataFrame(rows)
    table.to_csv(OUT/"station_trends.csv",index=False)
    print(table.round(3).to_string(index=False))
    for plant,g in panel.groupby("plant"):
        w=g.pivot(index="year",columns="land",values="onset_midpoint").dropna()
        print(plant,"paired onset",len(w),list(w.index),list(w.island-w.mainland))
    assert len(panel)==56 and not panel.duplicated(["plant","land","year"]).any()


if __name__=="__main__":
    run(Path(sys.argv[1]))
