"""Plot the four dated cuts using cloud-screened Sentinel-2 observations."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from alan_wet_shore import repo
out=repo/"results/all_species_forest/wet_annual"
v=pd.read_csv(out/"postcut_sentinel.csv.gz").dropna(subset=["ndvi"])
v=v.groupby(["y","since_cut"],as_index=False).ndvi.median()
fig,axes=plt.subplots(2,2,figsize=(8,5),sharex=True,sharey=True,layout="constrained")
for ax,(year,data) in zip(axes.flat,v.groupby("y")):
    ax.axvspan(0,14,color="#f8ebdb");ax.axvline(0,color="#9a6037",lw=.9)
    ax.scatter(data.since_cut,data.ndvi,color="#166358",s=19)
    ax.set_title(str(year),loc="left");ax.set_ylim(.05,1);ax.set_xlim(-30,44)
    ax.grid(axis="y",alpha=.18)
fig.supxlabel("Days from satellite-detected first cut")
fig.supylabel("Sentinel-2 NDVI (10 m)")
fig.suptitle("Vegetation before and after mowing at the Somerset wetland site",fontsize=13)
fig.savefig(out/"postcut.png",dpi=165)
