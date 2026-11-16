"""Report direction for the species with >=20 coordinates, and audit field onset for directional leads."""
from pathlib import Path
from zipfile import ZipFile
import sys
import duckdb
import numpy as np
import pandas as pd
import rasterio

cache, rasters, lights, out = map(Path, sys.argv[1:])
db = duckdb.connect(str(cache / "global_gbif.duckdb"))
summary = pd.read_csv(out / "species_screen.csv")
rows = []
for r in summary.query("coordinates >= 20 and blocks >= 10").itertuples():
    g = db.execute("select * from plant_coordinate_change where sp=?", [r.sp]).df()
    x = g[["delta", "light0", "rich_change", "change"]]
    a = (x - x.groupby(g.block).transform("mean")).to_numpy()
    resid = a[:,0] - a[:,1:3] @ np.linalg.lstsq(a[:,1:3],a[:,0],rcond=None)[0]
    info = pd.Series(resid**2).groupby(g.block).sum()
    bright = g.delta > 0
    rows.append(dict(species=r.species,coordinates=len(g),island=r.island_coordinates,
                     mainland=r.mainland_coordinates,brightened=int(bright.sum()),
                     mean_change_bright=g.change[bright].mean(),mean_change_other=g.change[~bright].mean(),
                     descriptive_adjusted_per_10dn=10*(resid@a[:,3])/(resid@resid),
                     informative_blocks=int((info>max(1e-10,info.sum()*1e-9)).sum())))
broad=pd.DataFrame(rows)
broad.to_csv(out/"high_coverage_descriptive.csv",index=False)
targets=broad.query("mean_change_bright<0 and descriptive_adjusted_per_10dn<0").species
with ZipFile(cache/"chronicle-of-nature-calendar.v1.0.5.zip") as archive:
    def table(name):
        member=next(m for m in archive.namelist() if m.endswith("/"+name+".csv"))
        return pd.read_csv(archive.open(member),low_memory=False)
    phen=table("phenology")
    stations=table("studysites")
p=phen[phen.taxon.isin(targets)&phen.eventtype.eq("onset of blooming")&
       phen.quality.eq("OK")&phen.year.between(1992,2014)].merge(stations,on="studysite",validate="m:1")
p["year"]=p.year.astype(int)
p=p.rename(columns={"taxon":"species","dayofyear":"flower"})
pieces=[]
for year,g in p.groupby("year"):
    g=g.copy()
    for col,path in [("sos",rasters/f"GIMMS3g_v1_SOS_{year}_Method_4_30_90.tif"),
                     ("dn",next(lights.glob(f"Harmonized_DN_NTL_{year}_*.tif")))]:
        with rasterio.open(path) as s:
            g[col]=[float(v[0]) if not np.ma.getmaskarray(v)[0] else np.nan for v in
                    s.sample(g[["longitude","latitude"]].values,masked=True)]
    pieces.append(g)
panel=pd.concat(pieces).query("1<=sos<=240 and 0<=dn<=63").copy()
panel["gap"]=panel.flower-panel.sos
panel["absolute_gap"]=panel.gap.abs()
assert not panel.duplicated(["studysite","species","year"]).any()
xy=db.sql("""select distinct species,longitude,latitude from plant_coordinate_change c join names n
             on c.sp=n.specieskey join site_exposure s on s.unit_id=c.xy""").df()
panel=panel.merge(xy.assign(exact_gbif_coordinate=True),how="left",validate="m:1")
panel["exact_gbif_coordinate"]=panel.exact_gbif_coordinate.eq(True)
db.register("panel",panel)
db.execute("create or replace table exploratory_flowering as select * from panel")
s=panel.groupby(["species","studysite"]).agg(years=("year","nunique"),start=("year","min"),end=("year","max"),
     dn_min=("dn","min"),dn_max=("dn","max"),mean_gap=("gap","mean"),min_gap=("absolute_gap","min"),
     max_gap=("absolute_gap","max"),exact_gbif=("exact_gbif_coordinate","any"))
s.to_csv(out/"exploratory_field_coverage.csv")
print(broad.to_string(index=False));print(s.to_string())
