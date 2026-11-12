"""Compare same-species region-year flowering series with four sampled light/land arms."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

GROUP = ["plant", "region5", "year"]
ARM = ["land_class", "light"]
SITE = ["source","station"]
UNIT = GROUP + ARM + SITE + ["individual"]
ARMS = [("island","bright"), ("island","dark"), ("mainland","bright"), ("mainland","dark")]


def four(frame):
    return (frame.groupby(GROUP + ARM).size().unstack(ARM)
            .reindex(columns=pd.MultiIndex.from_tuples(ARMS), fill_value=np.nan).dropna())


def read_source(name, visits, exposure):
    v = pd.read_parquet(visits,columns=["year","latitude","longitude","plant","station","individual","date","status"])
    e = pd.read_parquet(exposure).dropna(subset=["light"])
    v = v.merge(e, on=["year","latitude","longitude"],validate="m:1")
    v["source"] = name
    v["station"] = v.station.astype(str)
    v["marked"] = v.individual.notna()
    v["individual"] = v.individual.fillna("unmarked population at reported location").astype(str)
    v["day"] = pd.to_datetime(v.date).dt.dayofyear
    assert set(v.status.unique()) <= {0,1}
    return v


def extract_units(v):
    daily = v.groupby(UNIT + ["day"]).status.agg(["min","max"]).reset_index()
    consistent = daily[daily["min"] == daily["max"]].sort_values("day")
    summary = []
    for ix, g in consistent.groupby(UNIT):
        yes = g.loc[g["max"].eq(1), "day"].to_numpy()
        no = g.loc[g["max"].eq(0), "day"].to_numpy()
        if not len(yes) or not len(no):
            continue
        row = dict(zip(UNIT, ix))
        y0, y1 = yes[0], yes[-1]
        before, after = no[no < y0], no[no > y1]
        if len(before) and len(after) and not ((no > y0) & (no < y1)).any():
            row.update(lo=before[-1]+1, first=y0, last=y1, hi=after[0]-1)
        summary.append(row)
    return pd.DataFrame(summary, columns=UNIT + ["lo","first","last","hi"])


def outcomes(complete):
    rows = []
    for ix,g in complete.groupby(GROUP + ARM + SITE):
        row = dict(zip(GROUP + ARM + SITE,ix))
        row["duration_low"] = (g["last"]-g["first"]).mean()
        row["duration_high"] = (g.hi-g.lo).mean()
        if g.individual.nunique() >= 2:
            days = np.arange(int(g.lo.min()), int(g.hi.max())+1)
            certain = ((g["first"].to_numpy()[:,None] <= days) &
                       (days <= g["last"].to_numpy()[:,None])).sum(axis=0)
            possible = ((g.lo.to_numpy()[:,None] <= days) & (days <= g.hi.to_numpy()[:,None])).sum(axis=0)
            eligible = days[possible >= certain.max()]
            row["peak_low"],row["peak_high"] = min(eligible),max(eligible)
        rows.append(row)
    return pd.DataFrame(rows, columns=GROUP + ARM + SITE + ["duration_low","duration_high","peak_low","peak_high"])


def contrast(stations, metric):
    cols = [metric + "_low", metric + "_high"]
    rows=[]
    for ix,g in stations.dropna(subset=cols).groupby(GROUP):
        by = g.groupby(ARM)[cols].mean()
        if not set(ARMS) <= set(by.index):
            continue
        plus = by.loc[[ARMS[0],ARMS[3]]]
        minus = by.loc[[ARMS[1],ARMS[2]]]
        row = dict(zip(GROUP,ix))
        row.update(metric=metric, low_days=plus.iloc[:,0].sum()-minus.iloc[:,1].sum(),
                   high_days=plus.iloc[:,1].sum()-minus.iloc[:,0].sum(),
                   min_stations=int(g.groupby(ARM).size().min()))
        rows.append(row)
    return rows


def summarize(v, out):
    out.mkdir(exist_ok=True,parents=True)
    coverage_by_source = v.groupby("source").agg(records=("status","size"), species=("plant","nunique"), regions=("region5","nunique"))
    coverage_by_source.to_csv(out/"source_coverage.csv")
    dated = v.groupby(UNIT).day.nunique().rename('visit_days').reset_index()
    repeated = dated[dated.visit_days.ge(2)]
    units = extract_units(v)
    arm_status = v.groupby(GROUP+ARM).status.nunique().rename("state_count").reset_index()
    arm_both = arm_status[arm_status.state_count.eq(2)]
    complete = units.dropna(subset=["lo","hi"]).copy()
    station_outcomes = outcomes(complete)
    coverage = []
    stages = [("dated_binary",v),("at_least_two_dates",repeated),("both_states_across_locations",arm_both),("both_states_in_repeated_unit",units),
                        ("bounded_flowering_bout",complete),("population_peak_window",station_outcomes.dropna(subset=["peak_low"]))]
    for stage,group in stages:
        counts = four(group)
        for ix,row in counts.iterrows():
            coverage.append(dict(zip(GROUP,ix)) | dict(stage=stage, **{"_".join(k):int(row[k]) for k in ARMS}))
    cov=pd.DataFrame(coverage,columns=GROUP+["stage"]+["_".join(k) for k in ARMS])
    cov.to_csv(out/"four_arm_support.csv", index=False)

    detail = v.merge(cov[cov.stage.eq("dated_binary")][GROUP],on=GROUP).groupby(GROUP+ARM).agg(days=("day","nunique"),visits=("day","size"),open_visits=("status","sum")).reset_index()
    detail["visits_without_open_flowers"] = detail.visits-detail.open_visits
    detail=detail.merge(units.groupby(GROUP+ARM).size().rename("units_with_both_states"),on=GROUP+ARM,how="left")
    detail["units_with_both_states"]=detail.units_with_both_states.fillna(0).astype(int)
    detail.to_csv(out/"candidate_arm_counts.csv",index=False)
    estimates = pd.DataFrame(contrast(station_outcomes,"duration") + contrast(station_outcomes,"peak"),
                             columns=GROUP + ["metric","low_days","high_days","min_stations"])
    estimates.to_csv(out/"response_intervals.csv",index=False)
    print(cov.groupby("stage").agg(strata=("plant","size"),species=("plant","nunique"),regions=("region5","nunique")))
    print(estimates.to_string(index=False))
    totals = {stage:dict(strata=len(g),species=g.plant.nunique(),regions=g.region5.nunique())
              for stage,_ in stages for g in [cov[cov.stage.eq(stage)]]}
    totals["any_arm"] = dict(records=len(v),species=v.plant.nunique(),regions=v.region5.nunique(),region_years=len(v[["region5","year"]].drop_duplicates()))
    (out/"regional_counts.json").write_text(json.dumps(totals,indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",action="append",nargs=3,metavar=("NAME","VISITS","EXPOSURE"),required=True)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--frozen-sites",type=Path,help="Private 10,698-row table with frozen key and land_id")
    a=p.parse_args()
    joined = pd.concat([read_source(*r) for r in a.source],ignore_index=True)
    summarize(joined,a.out)
    if a.frozen_sites:
        frozen = pd.read_parquet(a.frozen_sites)
        assert len(frozen) == frozen.key.nunique() == 10698
        matching = joined.merge(frozen[["key","land_id","original_light"]],on=["key","land_id"],validate="m:1")
        matching["light"] = matching.original_light
        summarize(matching,a.out/"exact_cells")
