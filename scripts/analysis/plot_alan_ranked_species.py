"""Draw selected species and time contributions from the public summary tables."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

folder = Path(__file__).resolve().parents[2] / "results/focal_species"
species = pd.read_csv(folder / "ranked_contributions.csv").set_index("plant")
annual = pd.read_csv(folder / "annual_contributions.csv")
annual.loc[annual.island_site_years + annual.mainland_site_years == 0, "gap_contribution_pp"] = np.nan
annual = annual.pivot(index="plant", columns="year", values="gap_contribution_pp").loc[species.index]
y = np.arange(len(species))
fig, (left, right) = plt.subplots(1, 2, figsize=(10.8, 4.4), layout="constrained",
                                 gridspec_kw={"width_ratios": [1.12, 1.02]})
left.axvline(0, color="#9a9a9a", lw=0.8)
left.hlines(y, species.island_pp, species.mainland_pp, color="#b0b5b8", linewidth=2)
left.scatter(species.island_pp, y, label="Island", c="#006963", zorder=3, s=49)
left.scatter(species.mainland_pp, y, label="Mainland", c="#bd5c22", marker="D", zorder=3, s=38)
left.set_yticks(y, species.index, fontstyle="italic")
left.set_ylim(len(y) - 0.5, -0.5)
left.set_xlim(-0.21, 0.36)
left.set_xlabel("Percentage points per unit log(1 + VIIRS)")
left.set_title("Contribution to each land-class slope", loc="left", fontsize=12)
left.legend(loc="upper left", bbox_to_anchor=(0, -0.15), frameon=False, ncol=2)
left.spines[["top", "right", "left"]].set_visible(False)
right.set_title("Year's contribution to the pooled gap", loc="left", fontsize=12)
scale = max(abs(np.nanmin(annual)), abs(np.nanmax(annual)))
palette = plt.get_cmap("BrBG").copy()
palette.set_bad("#e4e6e7")
right.imshow(annual, aspect="auto", cmap=palette, norm=TwoSlopeNorm(0, -scale, scale))
right.set_xticks(range(len(annual.columns)), annual.columns)
right.set_yticks([])
right.set_xlabel("Grey cell = no record on either side")
for row in range(len(annual)):
    for col in range(len(annual.columns)):
        value = annual.iloc[row, col]
        label = "—" if np.isnan(value) else f"{value:+.3f}"
        right.text(col, row, label, va="center", ha="center", fontsize=8,
                   color="white" if abs(value) > scale * .70 else "#18201e")
for ax in (left, right):
    ax.axhline(2.5, color="#989d9a", lw=.8)
    ax.tick_params(length=0, pad=5)
fig.savefig(folder / "species_and_years.png", dpi=180, bbox_inches="tight")
plt.close(fig)
