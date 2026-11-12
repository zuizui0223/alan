"""Describe negative species contributions in the frozen model without fitting a new model."""
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[2] / "results"
data = pd.read_csv(root / "species_interpretation/shared_544.csv")
additions = pd.read_csv(root / "mechanism_exemplars/transport_updates.csv")
assert len(data) == data.plant.nunique() == 544
np.testing.assert_allclose(data.island_pp - data.mainland_pp, data.gap_pp, atol=1e-12)
negative = data.loc[data.gap_pp.lt(0)].copy()
assert len(negative) == 302 and not negative[["island_pp", "mainland_pp"]].eq(0).any().any()
negative["slope_pattern"] = np.select(
    [negative.island_pp.gt(0), negative.mainland_pp.gt(0)],
    ["both_positive_mainland_steeper", "island_negative_mainland_positive"],
    default="both_negative_island_steeper",
)
negative["transport_documented"] = negative.transport_review_status.eq("documented") | negative.plant.isin(additions.plant)
negative["rank_negative_gap"] = negative.gap_pp.rank(method="min").astype(int)
columns = ["rank_negative_gap", "plant", "slope_pattern", "island_pp", "mainland_pp", "gap_pp",
           "island_site_years", "mainland_site_years", "region_years_with_species_on_both_sides",
           "v011", "transport_documented", "night_exclusion_studies", "day_only_relative_loss",
           "glopl_cases", "glopl_median_log_response_ratio", "sc_source_class"]
negative = negative[columns].sort_values("rank_negative_gap")
negative.to_csv(root / "species_interpretation/negative_302.csv", index=False)
for pattern, group in negative.groupby("slope_pattern"):
    print(pattern, len(group), round(group.gap_pp.sum(),3), group.transport_documented.sum(),
          group.glopl_cases.gt(0).sum(), group.day_only_relative_loss.notna().sum(),
          group.region_years_with_species_on_both_sides.gt(0).sum(), group.v011.sum(),
          group.island_site_years.eq(1).sum())
for col in ["glopl_median_log_response_ratio","day_only_relative_loss"]:
    matched = negative.dropna(subset=[col])
    print(col, len(matched), matched[col].rank().corr(matched.gap_pp.rank()))
print(negative.dropna(subset=["day_only_relative_loss"]).to_string(index=False))
