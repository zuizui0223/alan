"""Replay denominator and geographic checks against screened adult-moth predation sources."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def verified_excel(path, digest, **kwargs):
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"Source workbook changed: {path.name}")
    return pd.read_excel(path, **kwargs)


def run(frame, japan, europe, exposure, out):
    out.mkdir(parents=True, exist_ok=True)
    sites = pd.read_parquet(frame).drop_duplicates("key")
    pairs = sites[["region5","year"]].drop_duplicates().sort_values(["region5","year"])
    assert (len(sites),len(pairs),pairs.region5.nunique()) == (10698,195,78)
    assert set(pairs.year) == {2019,2021,2023,2024,2025}
    jp = verified_excel(japan,"1bfd572eaf3be845c006f93c7f0670bf897871fe02fa7797be14fa2ed215b6b1",sheet_name=None)
    trial = jp["Predator attack"]
    moths = trial[trial.Visitor_group.eq("Moths")].copy()
    assert set(moths.Date.dt.year) == {2020}
    counts = moths.assign(attacked=moths.Attack.eq("Yes"),caught=moths.Predation.eq("Yes")).groupby("Period").agg(
        encounters=("Attack","size"),attacks=("attacked","sum"),captures=("caught","sum"))
    counts.to_csv(out/"japan_2020_mantis_encounters.csv")
    assert tuple(counts.loc["Night"]) == (39,16,5)
    # Named Honshu localities in the 2019 event sheet straddle the 135 E meridian.
    assert not pairs.region5.isin(["62_24","63_24"]).any()
    events = jp["Predation events"]
    n_2019 = len(events[events.Date.dt.year.eq(2019) & events.Prey_order.eq("Lepidoptera")])
    eu = verified_excel(europe,"4294e27d94a5f748f4e3c6723fab0154baf302525dc1e1081545a340d2f71195")
    assert eu.Day.dt.year.eq(2025).all() and eu.Catch_automated.notna().all()
    assert not eu.duplicated(["country","Trap code","Day"]).any()
    totals = eu.groupby(["country","habitat"]).agg(traps=("Trap code","nunique"),
            trap_days=("Day","size"),moths_caught_by_trap=("Catch_automated","sum")).reset_index()
    # Event counts are explicit in the Results of Bourougaaoui et al., version 1.
    labels = pd.DataFrame([
        ("Por","forest",38.835777,-8.849344,5,"forest_total"),
        ("Por","urban",38.713961,-9.188228,None,"event_count_unreported"),
        ("It","forest",46.630154,10.805001,3,"event_trap_subset_unidentified")],
        columns=["country","habitat","latitude","longitude","reported_fed_on","coverage"])
    centroids = pd.read_parquet(exposure)
    rows = labels.merge(centroids,on=["latitude","longitude"],validate="1:1").merge(totals,on=["country","habitat"],validate="1:1")
    assert len(rows)==3 and rows.valid_months.eq(4).all() and rows.land_class.eq("mainland").all()
    assert len(rows.merge(pairs,on=["region5","year"])) == 3
    rows["reported_fraction_all_trap_catches"] = (rows.reported_fed_on/rows.moths_caught_by_trap).where(rows.coverage.eq("forest_total"))
    rows.drop(columns=["latitude","longitude","key","land_id","unit_id","area_km2"]).to_csv(out/"trap_feeding.csv",index=False)
    # Counts refer to source groups that contain evidence, never zero predator attacks.
    by = rows[rows.reported_fed_on.notna()].groupby(["region5","year"]).size().rename("trap_groups_with_reported_feeding")
    pairs = pairs.merge(by,on=["region5","year"],how="left")
    pairs["trap_groups_with_reported_feeding"] = pairs.trap_groups_with_reported_feeding.fillna(0).astype(int)
    # Detailed candidate series from the source audit; retained rows need a directly scored
    # animal outcome and an exposure denominator. Study limitations are recorded alongside.
    source = pd.read_csv(out/"field_source_screen.csv",dtype={"years":str})
    natural = source[source.prey_context.isin(["free_adults","field_releases"]) & source.outcome.eq("direct_outcomes")].copy()
    natural["year"] = natural.years.str.split("|",regex=False)
    natural["region5"] = natural.regions.str.split("|",regex=False)
    natural = natural.explode("year").explode("region5")
    natural["year"] = natural.year.astype(int)
    matches = pairs[["region5","year"]].merge(natural,on=["region5","year"])
    if matches[["land_class","light"]].isna().any().any():
        raise ValueError("Classify light and land at newly matched field stations before counting four-arm support.")
    arms = matches.dropna(subset=["land_class","light"]).drop_duplicates(["region5","year","land_class","light"])
    n_arms = arms.groupby(["region5","year"]).size().rename("supported_direct_rate_arms")
    pairs = pairs.merge(n_arms,on=["region5","year"],how="left")
    pairs["supported_direct_rate_arms"] = pairs.supported_direct_rate_arms.fillna(0).astype(int)
    pairs["free_adult_moth_four_arm_rate_available"] = pairs.supported_direct_rate_arms.eq(4)
    pairs.to_csv(out/"region_year_coverage.csv",index=False)
    result = dict(regions_screened=78,region_years_screened=195,reported_trap_feeding_regions=int(rows[rows.reported_fed_on.notna()].region5.nunique()),
                  free_moth_four_arm_regions=int(pairs[pairs.free_adult_moth_four_arm_rate_available].region5.nunique()),
                  jp_2019_lepidoptera_events=n_2019, estimable_island_minus_mainland_light_slope=None)
    (out/"summary.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ["frame","japan","europe","exposure","out"]:
        p.add_argument(arg,type=Path)
    a=p.parse_args()
    run(a.frame,a.japan,a.europe,a.exposure,a.out)
