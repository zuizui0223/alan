"""Verify the species ledger and display out-of-block light importance in the two splits."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import duckdb
import geopandas as gpd


def report(work,results):
    root=results/"all_species_forest"
    old=set(pd.read_csv(results/"long_term_land_contrasts/environment_forest/species.csv").sp)
    previous=pd.read_csv(results/"long_term_land_contrasts/environment_forest/shortlist.csv")[["sp","species"]]
    db=duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    xy=db.sql("select distinct unit_id,longitude,latitude from geo").df()
    world=gpd.read_file(work/"gbif_annual/countries.zip")[["CONTINENT","ADMIN","geometry"]]
    points=gpd.GeoDataFrame(xy,geometry=gpd.points_from_xy(xy.longitude,xy.latitude),crs=4326)
    geo=gpd.sjoin(points,world,how="left",predicate="within")[["unit_id","CONTINENT","ADMIN"]].fillna("unassigned")
    assert geo.unit_id.is_unique
    rows=[]
    for name,folder in [("1992–2024",root),("1992–2014, green-up",root/"greenup")]:
        ledger=pd.read_csv(folder/"all_species.csv.gz")
        repeated=pd.read_csv(folder/"repartition.csv")
        undetermined=repeated.light_response.isna() | repeated.light_response.eq(0)
        if (undetermined & repeated.direction.ne("undetermined")).any():
            repeated.loc[undetermined,"direction"]="undetermined"
            repeated.to_csv(folder/"repartition.csv",index=False)
        candidates=pd.read_csv(folder/"candidates.csv")
        assert len(ledger)==148786 and ledger.sp.is_unique
        fit=ledger[ledger.status.eq("evaluated")]
        selected=fit.top.eq("light") & fit.skill.gt(0) & fit.gain_light.gt(0)
        assert set(fit.loc[selected,"sp"])==set(repeated.sp)
        assert set(candidates.sp)==set(repeated.loc[repeated.stable_light,"sp"])
        raw=pd.read_parquet(work/"gbif_annual"/("all_species_rf_sites.parquet" if folder==root else "green_rf_sites.parquet"))
        places=raw[["unit_id","island"]].drop_duplicates().merge(geo,on="unit_id",validate="1:1")
        places.groupby(["CONTINENT","ADMIN","island"]).size().rename("coordinates").to_csv(folder/"geography.csv")
        assert set(raw.sp)==set(ledger.loc[ledger.status.ne("no_qualifying_long_series"),"sp"])
        matched=candidates[candidates.net_agrees].copy()
        rows.append(dict(window=name,roster=len(ledger),with_series=raw.sp.nunique(),sites=raw.unit_id.nunique(),
                         spatial_evaluated=len(fit),first_light=len(repeated),stable=len(candidates),
                         stable_positive=candidates.direction.eq("positive").sum(),stable_negative=candidates.direction.eq("negative").sum(),
                         observed_increase=matched.direction.eq("positive").sum(),observed_decrease=matched.direction.eq("negative").sum(),
                         stable_outside_previous_screen=(~candidates.sp.isin(old)).sum(),
                         agreeing_outside_previous_screen=(~matched.sp.isin(old)).sum()))
        candidates["smaller_gain"]=candidates[["gain_light","gain_light_repeat"]].min(axis=1)
        candidates["in_previous_427"]=candidates.sp.isin(old)
        ranked=candidates.sort_values(["net_agrees","smaller_gain"],ascending=False)
        ranked.to_csv(folder/"ranked.csv",index=False)
        ranked[ranked.net_agrees].to_csv(folder/"observed_and_predicted.csv",index=False)
    pd.DataFrame(rows).to_csv(root/"coverage.csv",index=False)
    ranked=pd.read_csv(root/"ranked.csv")
    ledger=pd.read_csv(root/"all_species.csv.gz")
    previous.merge(ledger[["sp","status","top","gain_light","skill","light_response"]],on="sp",
                   validate="1:1").merge(ranked[["sp","direction","net_agrees"]],on="sp",how="left",
                   validate="1:1").to_csv(root/"previous_four.csv",index=False)
    shortlist=pd.concat([g.nlargest(6,"smaller_gain") for _,g in ranked[ranked.net_agrees].groupby("direction")])
    if shortlist.empty:return
    shortlist=shortlist.reset_index(drop=True)
    fig,ax=plt.subplots(figsize=(7,max(3,.34*len(shortlist)+1.6)))
    for i,r in shortlist.iterrows():
        a,b=100*r.gain_light,100*r.gain_light_repeat
        color="#277285" if r.direction=="negative" else "#a94b21"
        ax.plot([a,b],[i,i],color=color,lw=2)
        ax.scatter([a,b],[i,i],c=["white",color],edgecolors=color,s=35)
    ax.set_yticks(range(len(shortlist)),[f"{r.species}  ({'decrease' if r.direction=='negative' else 'increase'})" for r in shortlist.itertuples()])
    ax.invert_yaxis()
    ax.set_xlabel("Prediction error added by omitting light (% of baseline MSE)")
    ax.set_xlim(left=0)
    fig.text(.02,.01,"Open / filled circles: initial / second spatial split. Up to six per direction.",fontsize=8)
    fig.tight_layout(rect=(0,.045,1,1))
    fig.savefig(root/"light_candidates.png",dpi=160)
    plt.close(fig)


if __name__=="__main__":
    report(*map(Path,sys.argv[1:3]))
