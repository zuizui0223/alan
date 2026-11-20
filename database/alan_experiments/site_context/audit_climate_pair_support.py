"""Check spatial resolution before using climate as a within-pair covariate."""
from pathlib import Path
import hashlib
import json
import math
import pandas as pd

OUT=Path(__file__).resolve().parent
SOURCE=OUT/"annual_climate.csv"

def distance(a,b,c,d):
    a,b,c,d=map(math.radians,[a,b,c,d])
    return 12742*math.asin(min(1,math.sqrt(math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2)))

def main():
    climate=pd.read_csv(SOURCE)
    g=climate[climate.study_id.eq("Giavi2021")].copy()
    g["pair"]=g.site_id.str.extract(r"pair(\d+)")[0].astype(int)
    g["arm"]=g.site_id.str.extract(r"_(dark|illuminated)$")[0]
    assert len(g)==12 and set(g.pair)==set(range(1,7))
    assert not g.duplicated(["pair","arm","year"]).any()
    rows=[]
    for (pair,year),s in g.groupby(["pair","year"]):
        assert set(s.arm)=={"dark","illuminated"}
        lit=s[s.arm.eq("illuminated")].iloc[0];dark=s[s.arm.eq("dark")].iloc[0]
        same=(lit.grid_lat,lit.grid_lon)==(dark.grid_lat,dark.grid_lon)
        dt=float(lit.temperature_mean_C-dark.temperature_mean_C)
        dp=float(lit.precipitation_total_mm-dark.precipitation_total_mm)
        if same:
            assert abs(dt)<1e-12 and abs(dp)<1e-10
        rows.append(dict(pair=int(pair),year=int(year),lit_site_id=lit.site_id,dark_site_id=dark.site_id,
            lit_grid_lat=lit.grid_lat,lit_grid_lon=lit.grid_lon,dark_grid_lat=dark.grid_lat,dark_grid_lon=dark.grid_lon,
            requested_site_separation_km=distance(lit.requested_lat,lit.requested_lon,dark.requested_lat,dark.requested_lon),
            grid_centre_separation_km=distance(lit.grid_lat,lit.grid_lon,dark.grid_lat,dark.grid_lon),same_grid=same,
            lit_minus_dark_temperature_C=dt,lit_minus_dark_precipitation_mm=dp,
            interpretation="shared_coarse_context_not_measured_equality" if same else "grid_boundary_contrast_not_verified_plot_climate"))
    out=pd.DataFrame(rows);out.to_csv(OUT/"climate_pair_support.csv",index=False)
    audit=dict(pairs=len(out),same_grid_pairs=int(out.same_grid.sum()),boundary_pairs=out.loc[~out.same_grid,"pair"].tolist(),
               annual_source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
               conclusion="No within-pair climate adjustment adopted. Neither equal-grid equality nor boundary difference proves actual plot climate. All input sites retained; no reanalysis or imputation.")
    (OUT/"climate_pair_audit.json").write_text(json.dumps(audit,indent=2)+"\n",encoding="utf-8")
    print(out.to_string(index=False))

if __name__=="__main__":
    main()
