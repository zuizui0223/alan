"""Show annual observations with missing years left open."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from alan_long_phenology import OUT

panel=pd.read_csv(OUT/"station_year.csv")
fig,axes=plt.subplots(3,2,figsize=(9,7),sharex=True,layout="constrained")
for col,(plant,group) in enumerate(panel.groupby("plant")):
    axes[0,col].set_title(plant,fontstyle="italic")
    for row,(field,label) in enumerate([("onset_midpoint","Flower onset (day of year)"),
                                        ("ndvi_on_open_dates","NDVI near observed flowering"),
                                        ("viirs","Night radiance (nW cm⁻² sr⁻¹)")]):
        ax=axes[row,col]
        for side,color,marker in [("island","#276f61","o"),("mainland","#a14d2c","s")]:
            g=group[group.land.eq(side)].sort_values("year")
            ax.plot(g.year,g[field],color=color,marker=marker,markersize=4,label=side,linewidth=1)
        if col==0: ax.set_ylabel(label)
        if row==2: ax.set_xlabel("Year")
        ax.set_xticks([2012,2016,2020,2024])
        ax.grid(axis="y",alpha=.2)
axes[0,0].legend(frameon=False)
fig.savefig(OUT/"annual_series.png",dpi=180)
