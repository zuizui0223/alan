"""Replay the PR 85 species slopes, then attach independently documented autonomy."""
from pathlib import Path
import sys

import pandas as pd

source = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(source))
from scripts.analysis import alan_diel_repro_species_ghsl_source_gate_20261116 as pr

output = Path(__file__).resolve().parents[2] / "results/joint_night_autonomy"
sites = pr.load_sites(next((source / "pcnl/joined").glob("*PCNLPrior*.csv.gz")))
flora = pr.original_accepted_species(source / "groups", sites, pr.load_labels(Path(sys.argv[2])))
traits = pd.read_csv(output / "species.csv").set_index("species")
matched = traits.loc[traits.index.isin(flora.species.unique())].copy()
assert len(matched) == 52 and matched.night_stronger.sum() == 27
population, _ = pr.source_annual_population(sites, source / "annual")
photo, _ = pr.within_year_island_photo_design(population)
climate, _ = pr.matching_climate_gradient(photo, source / "climate")
built, _ = pr.authentic_ghsl_source(source / "ghsl/ghsl2015_2020", climate)
assert len(built) == 64944 and built.island.nunique() == 202
guild = matched.night_stronger.map({True: "night_higher", False: "day_higher"})
positive, strata, _ = pr.candidate_species_photo_support(built, flora, guild)
for name in matched.index:
    if name not in positive or name not in strata:
        matched.loc[name, "fit_status"] = "no informative island-year"
        continue
    result, reason = pr.fit_source_species(built, positive[name], strata[name])
    if reason:
        matched.loc[name, "fit_status"] = reason
        continue
    stats = result[0]
    matched.loc[name, "fit_status"] = "estimated"
    for col in ["n_islands", "n_photo_species_positive", "L_built2015_pp"]:
        matched.loc[name, col] = stats[col]
fitted = matched.dropna(subset=["L_built2015_pp"])
assert len(fitted) == 39 and (fitted.n_islands >= 10).sum() == 22
assert round(fitted.groupby("night_stronger").L_built2015_pp.mean().diff().iloc[-1], 3) == -1.646
matched.to_csv(output / "pr85_island_species.csv")
print(matched[matched.night_stronger & matched.autonomy.eq("absent")][["night_share", "n_islands", "L_built2015_pp", "fit_status"]])
