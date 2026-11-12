from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[2]
folder = root / "results/species_interpretation"
species = pd.read_csv(folder / "shared_544.csv").set_index("plant")
updates = pd.read_csv(root / "results/mechanism_exemplars/transport_updates.csv")
species.loc[updates.plant, "transport_review_status"] = "documented"
island, mainland = species.island_pp.lt(0), species.mainland_pp.lt(0)
groups = {"both_negative": island & mainland, "island_negative_only": island & ~mainland,
          "mainland_negative_only": ~island & mainland, "both_positive": ~island & ~mainland}
assert species.index.is_unique and len(species) == 544 and (species[["island_pp","mainland_pp"]] != 0).all().all()
rows = []
for group, mask in groups.items():
    cohort = species[mask]
    pl = cohort.glopl_median_log_response_ratio.dropna()
    night = cohort.day_only_relative_loss.dropna()
    rows.append(dict(group=group, species=len(cohort), moth_transport_documented=cohort.transport_review_status.eq("documented").sum(),
                     pl_measured=len(pl), pl_positive=(pl > 0).sum(), pl_median=pl.median(),
                     night_measured=len(night), night_loss_positive=(night > 0).sum(),
                     paired_region_year_support=cohort.region_years_with_species_on_both_sides.gt(0).sum(),
                     island_site_years_median=cohort.island_site_years.median(), mainland_site_years_median=cohort.mainland_site_years.median()))
table = pd.DataFrame(rows)
assert table.species.sum() == 544 and table.pl_measured.sum() == 64 and table.night_measured.sum() == 12
table.to_csv(folder / "slope_direction_traits.csv", index=False)
print(table.to_string(index=False))
