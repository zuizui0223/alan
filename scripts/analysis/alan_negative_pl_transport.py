"""Audit overlap of frozen composition contrasts, pollen transport and dated GloPL cases."""
from pathlib import Path
import re

import pandas as pd

root = Path(__file__).resolve().parents[2] / "results"
traits = pd.read_csv(root / "species_interpretation/shared_544.csv")
updates = pd.read_csv(root / "mechanism_exemplars/transport_updates.csv")
cases = pd.read_csv(root / "moth_pollen_evidence/pollen_limitation_cases.csv")
assert traits.plant.is_unique and len(traits) == 544
assert len(cases) == 261 and cases.plant.nunique() == 64
target = {2019, 2021, 2023, 2024, 2025}

def study_years(value):
    if pd.isna(value):
        return set()
    assert re.fullmatch(r"\d{4}(-\d{4})?", str(value)), value
    bounds = list(map(int, str(value).split("-")))
    return set(range(bounds[0], bounds[-1] + 1))

years = cases.year_of_study.map(study_years)
assert years.map(bool).sum() == 202 and not years.map(lambda x: bool(x & target)).any()
known = set().union(*years)
assert (min(known), max(known)) == (1977, 2015)
traits["documented"] = traits.transport_review_status.eq("documented") | traits.plant.isin(updates.plant)
eligible = traits[traits.island_pp.lt(0) & traits.gap_pp.lt(0) & traits.glopl_cases.gt(0)]
assert len(eligible) == 36
columns = ["plant", "island_pp", "mainland_pp", "gap_pp", "v011", "region_years_with_species_on_both_sides",
           "glopl_cases", "glopl_median_log_response_ratio", "day_only_relative_loss", "documented"]
candidates = eligible[eligible.documented].loc[:, columns].sort_values("gap_pp")
assert len(candidates) == 12
out = root / "mechanism_exemplars/negative_pl_transport.csv"
candidates.to_csv(out, index=False)
