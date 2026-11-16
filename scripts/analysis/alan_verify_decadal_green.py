"""Verify archived green-up, then compare both MODIS definitions on identical observations."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as sm

from alan_decadal_light import C, OUT
from alan_extended11 import coords
from alan_review_satellite_spring import greenup_date


def fit(data, outcome, light, unit):
    frame = data.copy()
    frame["x"] = frame[light] * (1 / np.log(1.1) if light == "loglight" else 0.1)
    frame["ix"] = frame.x * frame.land.eq("island")
    frame["y"] = frame[outcome]
    model = sm.ols(f"y ~ x + ix + C({unit}) + C(year)", data=frame, eval_env=-1)
    assert np.linalg.matrix_rank(model.exog) == model.exog.shape[1]
    fit_ = model.fit(cov_type="cluster", cov_kwds={"groups": frame.station}, use_t=True)
    rows = []
    for side, term in [("island", "x+ix"), ("mainland", "x"), ("difference", "ix")]:
        test = fit_.t_test(term + "=0")
        lo, hi = test.conf_int()[0]
        rows.append(dict(outcome=outcome, light=light, side=side, days=test.effect.item(),
                         lo=lo, hi=hi, rows=len(frame), stations=frame.station.nunique()))
    return rows


def run():
    baseline = pd.read_csv(OUT / "green_source_batches.csv").set_index("file").sha256
    records = []
    manifest = []
    for station, point in coords.iterrows():
        for first in [2001, 2009, 2017]:
            path = C / f"green_{station}_{first}.json"
            digest = None
            if path.exists():
                raw = path.read_bytes()
                digest = hashlib.sha256(raw).hexdigest()
                source = json.loads(raw)
                assert source["nrows"] == source["ncols"] == 1
                assert np.allclose([source["latitude"],source["longitude"]],
                                   [point.latitude,point.longitude],rtol=0,atol=1e-5)
                bands = {}
                for entry in source["subset"]:
                    year = int(entry["modis_date"][1:5])
                    assert first <= year < first + 8
                    bands.setdefault(year, {})[entry["band"]] = entry["data"][0]
                for year, values in bands.items():
                    date = greenup_date(values)
                    if pd.notna(date) and year == date.year:
                        records.append(dict(station=station,year=year,green_mcd=date.dayofyear))
            if path.name in baseline:
                assert digest == baseline[path.name]
            manifest.append(dict(file=path.name,sha256=digest))
    out = OUT / "verification"
    out.mkdir(exist_ok=True)
    pd.DataFrame(manifest).to_csv(out / "coverage.csv",index=False)
    annual = pd.read_csv(OUT / "site_year.csv").merge(pd.DataFrame(records),on=["station","year"],how="left",validate="1:1")
    original = annual.dropna(subset=["green_doy"])
    np.testing.assert_array_equal(original.green_doy,original.green_mcd)
    annual[["station","year","land","green_mcd","light_dn"]].to_csv(out/"annual.csv",index=False)
    valid = annual.dropna(subset=["green_mcd","light_dn"])
    valid = valid[valid.groupby("station").year.transform("size")>=2]
    results = fit(valid,"green_mcd","light_dn","station")
    eras = []
    for label, lo, hi in [("DMSP",2001,2013),("VIIRS_equivalent",2014,2024)]:
        era = valid[valid.year.between(lo,hi)]
        era = era[era.groupby("station").year.transform("size")>=2]
        eras.extend(dict(era=label,**row) for row in fit(era,"green_mcd","light_dn","station"))
    pd.DataFrame(eras).to_csv(out/"eras.csv",index=False)
    prior = pd.read_csv(OUT.parent/"long_phenology_11/panel.csv")
    same = prior.merge(annual[["station","year","green_mcd","light_dn"]],on=["station","year"],validate="m:1")
    same = same.dropna(subset=["green","green_mcd","onset","light_dn","loglight"])
    same = same[same.groupby("unit").year.transform("size")>=2].copy()
    same["gap_mcd"] = same.onset - same.green_mcd
    for light in ["loglight","light_dn"]:
        for outcome in ["gap","gap_mcd","onset","green","green_mcd"]:
            results.extend(fit(same,outcome,light,"unit"))
    same.to_csv(out/"matched.csv",index=False)
    table = pd.DataFrame(results)
    table.to_csv(out/"estimates.csv",index=False)
    print("batches",pd.DataFrame(manifest).sha256.count(),"matched",same.groupby("land").size().to_dict())
    print(table.loc[table.side.eq("island")].round(2).to_string(index=False))


if __name__ == "__main__":
    run()
