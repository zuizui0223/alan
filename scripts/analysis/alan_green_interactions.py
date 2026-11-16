"""Compare woody stocks and observed animal-plant timing on verified green-up years."""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as sm

ROOT = Path(__file__).resolve().parents[2] / "results/decadal_phenology"
OUT = ROOT / "green_interactions"
OUT.mkdir(exist_ok=True)
annual = pd.read_csv(ROOT/"verification/annual.csv").dropna(subset=["green_mcd","light_dn"])
key = ["station","year","land"]
sites = pd.read_csv(ROOT/"site_year.csv")[key+["pixel"]]
env = annual.merge(sites,on=key,validate="1:1")
wood = pd.read_csv(ROOT/"biomass_pixel_year.csv")
env = env.merge(wood[["pixel","year","land","woody_tonnes"]],on=["pixel","year","land"],how="left",validate="m:1")
pairs = pd.read_csv(ROOT/"paired_dates.csv").merge(env,on=key,validate="m:1")
pairs["absolute_gap"] = pairs.animal_minus_plant_median_observed_day.abs()
pairs["unit"] = pairs.station.astype(str)+"|"+pairs.plant+"|"+pairs.animal
pairs.to_csv(OUT/"matched_pairs.csv",index=False)
bio = env.dropna(subset=["woody_tonnes"]).groupby(["pixel","year","land"],as_index=False).agg(
    woody_tonnes=("woody_tonnes","first"),green_mcd=("green_mcd","mean"),light_dn=("light_dn","mean"))
bio.to_csv(OUT/"matched_wood.csv",index=False)
estimates=[]
support=[]
def analyze(d,unit,cluster,outcome,label):
    for side,g in d.groupby("land"):
        a=g[g.groupby(unit).year.transform("nunique")>=2].copy()
        support.append(dict(group=label,side=side,matched=len(g),matched_stations=g[cluster].nunique(),
            repeated=len(a),clusters=a[cluster].nunique(),years=a.year.nunique()))
        if a[cluster].nunique()<3:continue
        a["light10"]=a.light_dn/10
        a["green10"]=a.green_mcd/10
        a["u"]=a[unit].astype(str)
        a["y"]=a[outcome]
        for predictors in ["light10","green10","light10 + green10"]:
            model=sm.ols("y ~ "+predictors+" + C(u) + C(year)",a)
            if np.linalg.matrix_rank(model.exog)<model.exog.shape[1] or model.df_resid<=0:continue
            fit=model.fit(cov_type="cluster",cov_kwds={"groups":a[cluster]},use_t=True)
            for term in predictors.split(" + "):
                lo,hi=fit.conf_int().loc[term]
                estimates.append(dict(group=label,side=side,outcome=outcome,predictors=predictors,
                    term=term,estimate=fit.params[term],low=lo,high=hi))
analyze(bio,"pixel","pixel","woody_tonnes","wood")
for guild,g in pairs.groupby("guild"):
    analyze(g,"unit","station","absolute_gap",guild)
    analyze(g,"unit","station","animal_minus_plant_median_observed_day",guild)
    analyze(g,"unit","station","green_mcd",guild+"_green")
pd.DataFrame(estimates).query("not (outcome == 'green_mcd' and predictors != 'light10')").to_csv(OUT/"estimates.csv",index=False)
pd.DataFrame(support).drop_duplicates().to_csv(OUT/"coverage.csv",index=False)
print(pd.DataFrame(support).drop_duplicates().to_string(index=False))
print(pairs.groupby(["land","guild"]).woody_tonnes.count())
