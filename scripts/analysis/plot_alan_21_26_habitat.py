"""Plot the 38 focal species-site-years and the dates on which the plants were observed."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator
import numpy as np
import pandas as pd

folder = Path(__file__).resolve().parents[2] / "results/local_21_26"
frame = pd.read_csv(folder / "focal_dated_sites.csv")
colors = {"Asclepias syriaca": "#87529d", "Cirsium vulgare": "#087765"}
markers = {"Asclepias syriaca": "o", "Cirsium vulgare": "^"}
fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), sharey=True, layout="constrained")
for col, side in enumerate(("island", "mainland")):
    for plant in colors:
        points = frame[frame.land_class.eq(side) & frame.plant.eq(plant)].copy()
        top, bottom = axes[:, col]
        top.scatter(points.photo_built_pct, np.log1p(points.viirs_median),
                    c=colors[plant], marker=markers[plant], s=46, alpha=.8)
        dates = points.assign(date=points.observation_dates.str.split("|", regex=False)).explode("date")
        day = pd.to_datetime("2000-" + dates.date.str[5:]).dt.dayofyear
        bottom.scatter(day, np.log1p(dates.viirs_median), c=colors[plant], marker=markers[plant],
                       s=34, alpha=.68)
    top.set_title(side.capitalize())
    top.set(xlim=(-5, 105), xlabel="Built area within 100 m of reported locations (%)")
    bottom.set(xlim=(59, 340), xlabel="Local observation date (five years pooled)")
    bottom.set_xticks([61, 122, 183, 245, 306], ["Mar", "May", "Jul", "Sep", "Nov"])
for ax in axes.flat:
    ticks = np.array([0, 1, 3, 10, 30, 100])
    ax.yaxis.set_major_locator(FixedLocator(np.log1p(ticks)))
    ax.yaxis.set_major_formatter(FixedFormatter([str(x) for x in ticks]))
    ax.grid(axis="y", alpha=.18)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
for ax in axes[:, 0]:
    ax.set_ylabel("VIIRS radiance (log1p scale)")
handles = [Line2D([], [], marker=markers[p], linestyle="none", label=p, color=colors[p]) for p in colors]
fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
fig.savefig(folder / "habitat_and_dates.png", dpi=180, bbox_inches="tight")
plt.close(fig)
