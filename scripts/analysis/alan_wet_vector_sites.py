"""Compare wet plants with insect versus wind vectors within the same surveyed coordinates."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t
import statsmodels.api as sm

repo=Path(__file__).resolve().parents[2]
out=repo/"results/all_species_forest/wet_reproduction"
cohort=pd.read_csv(out/"species.csv")
v=cohort["Pollination syndrome"].fillna("")
insect=v.str.contains("insect-pollination")
wind=v.str.contains("wind-pollination")
breed=cohort["Generative reproduction type"].fillna("")
strict=breed.isin(["allogamy","allogamy self-incompatibility"])
can_self=breed.str.contains("autogamy|mixed mating") & ~strict
source=pd.read_parquet(repo.parent/"gbif_annual/all_species_rf_sites.parquet")
controls=["temperature","precipitation","built","light0","temp_mean","rain_mean","built0","effort","log_days_trend","start_year","span","years","spring_shift","summer_shift","autumn_shift"]
env=["light",*controls]
keys=["unit_id","block","island"]
support=[]
for contrast,left,right,mask_left,mask_right in [
    ("insect_minus_wind","insect","wind",insect & ~wind,wind & ~insect),
    ("outcross_minus_self_capable_in_insect_group","outcross","self_capable",
     insect & ~wind & strict,insect & ~wind & can_self)]:
    selected=cohort[mask_left | mask_right].copy()
    selected["group"]=np.where(mask_left.loc[selected.index],left,right)
    raw=source.merge(selected[["sp","group"]],on="sp",validate="m:1")
    raw=raw.dropna(subset=["change",*env])
    sites=raw.groupby(keys)[env]
    valid=(sites.max()-sites.min()).lt(1e-7).all(axis=1)
    raw=raw.merge(valid[valid].index.to_frame(index=False),on=keys,validate="m:1")
    paired=raw.groupby([*keys,"group"]).change.mean().unstack().dropna()
    paired["rate_difference"]=paired[left]-paired[right]
    paired=paired.join(raw.groupby(keys)[env].first()).reset_index()
    
    for label,sites in [("all",paired),("island",paired[paired.island.eq(1)]),("mainland",paired[paired.island.eq(0)])]:
        sites=sites.groupby("block").filter(lambda g:len(g)>1 and g.light.max()-g.light.min()>1e-7).copy()
        members=raw[raw.unit_id.isin(sites.unit_id)].groupby("group").sp.nunique()
        record=dict(contrast=contrast,scope=label,coordinates=len(sites),blocks=sites.block.nunique(),
                    left_species=members.get(left,0),right_species=members.get(right,0))
        if record["blocks"]>1:
            cols=["rate_difference",*env]
            within=sites[cols]-sites.groupby("block")[cols].transform("mean")
            x=within[env].loc[:,within[env].std().gt(1e-9)].copy()
            x["light"]/=10
            model=sm.OLS(within.rate_difference,x).fit(cov_type="cluster",cov_kwds={"groups":sites.block,"use_correction":False})
            resid_df=len(sites)-sites.block.nunique()-np.linalg.matrix_rank(x)
            correction=np.sqrt(record["blocks"]/(record["blocks"]-1)*(len(sites)-1)/resid_df)
            half=t.ppf(.975,record["blocks"]-1)*model.bse["light"]*correction
            record.update(estimate=model.params["light"],low=model.params["light"]-half,high=model.params["light"]+half)
        support.append(record)
pd.DataFrame(support).to_csv(out/"within_site.csv",index=False)
