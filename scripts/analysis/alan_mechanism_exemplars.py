"""Locate experimental mechanism exemplars within the frozen photographic frame."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
# Deliberate case selections by measured process; no ranking by the fitted coefficient.
CASES = {
    "Silene latifolia": ("nursery_pollination_under_lighting", "10.1038/s41598-020-75471-1"),
    "Vaccinium angustifolium": ("night_service_and_daytime_redundancy", "10.1139/cjps-2021-0216"),
    "Clerodendrum trichotomum": ("visitor_frequency_vs_effective_pollination", "10.1371/journal.pone.0033803"),
    "Saponaria officinalis": ("measured_pollen_delivery_by_a_frequent_visitor", "10.1016/j.flora.2005.07.010"),
    "Datura wrightii": ("effective_hawkmoth_pollination", "10.1093/aob/mcp053"),
    "Veronica persica": ("seasonal_flower_availability_under_lighting", "10.1016/j.biocon.2025.111650"),
    "Heracleum maximum": ("adult_moth_prey_and_shared_flower_candidate", "10.1111/een.13399"),
    "Lonicera implexa": ("within_island_variation_in_floral_fit", "10.1093/aobpla/plz078"),
    "Plantago lanceolata": ("wind_pollination_and_moth_carriage", "10.6084/m9.figshare.c.4955030"),
}


def build(photo_file, output):
    ledger = pd.read_csv(ROOT / "results/species_interpretation/shared_544.csv").set_index("plant")
    assert len(ledger) == ledger.index.nunique() == 544
    photos = pd.read_parquet(photo_file)
    assert set(photos.plant) == set(ledger.index)
    sites = photos.drop_duplicates(["plant", "key"])
    sites_in_both = sites.groupby(["plant", "region5", "year"]).land_class.nunique()
    paired = sites_in_both[sites_in_both == 2].reset_index()
    out = []
    for plant, (mechanism, doi) in CASES.items():
        r = ledger.loc[plant]
        focal = photos[photos.plant.eq(plant)]
        pos = sites[sites.plant.eq(plant)].groupby("land_class").size()
        pairs = paired[paired.plant.eq(plant)]
        assert [pos.island, pos.mainland, len(pairs)] == [
            r.island_site_years, r.mainland_site_years, r.region_years_with_species_on_both_sides]
        months = focal.groupby(["region5", "year", "month"]).land_class.nunique()
        arms = focal.groupby(["region5", "year", "land_class", "light"]).size()
        four = arms.groupby(level=[0,1]).size().eq(4)
        out.append(dict(plant=plant, mechanism=mechanism, study_doi=doi, legacy_positive=r.legacy,
                        v011_positive=r.v011, island_site_years=pos.island, mainland_site_years=pos.mainland,
                        island_pp=r.island_pp, mainland_pp=r.mainland_pp, gap_pp=r.gap_pp,
                        paired_regions=pairs.region5.nunique(), paired_region_years=len(pairs),
                        paired_region_year_ids="|".join(pairs.region5 + ":" + pairs.year.astype(str)),
                        bilateral_region_year_months=int(months.eq(2).sum()),
                        region_years_with_all_four_photo_arms=int(four.sum())))
    frame = pd.DataFrame(out)
    np.testing.assert_allclose(frame.island_pp - frame.mainland_pp, frame.gap_pp)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("photos", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build(args.photos, args.output)
