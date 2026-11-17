"""Pool floral visitation with live crab spiders versus controls, then audit wet plant links."""
from pathlib import Path
import io
import subprocess

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import t

repo=Path(__file__).resolve().parents[2]
out=repo/"results/all_species_forest/spider_mechanism"
effects=pd.read_csv(out/"live_crab_spider.csv")
rows=[]
for rho in [0,0.5,1]:
    estimates=[]
    variances=[]
    for _,g in effects.groupby("source"):
        weights=1/g.variance
        weights/=weights.sum()
        estimates.append(weights@g.log_ratio)
        variances.append((1-rho)*((weights**2)@g.variance)+rho*(weights@np.sqrt(g.variance))**2)
    y,v=np.array(estimates),np.array(variances)
    df=len(y)-1
    def q(tau2):
        w=1/(v+tau2)
        return w@((y-np.average(y,weights=w))**2)
    tau2=0 if q(0)<=df else brentq(lambda x:q(x)-df,0,100)
    w=1/(v+tau2)
    mu=np.average(y,weights=w)
    half=t.ppf(.975,df)*np.sqrt(max(1,q(tau2)/df)/w.sum())
    rows.append(dict(rho=rho,publications=len(y),contrasts=len(effects),tau2=tau2,
                     ratio=np.exp(mu),low=np.exp(mu-half),high=np.exp(mu+half)))
pd.DataFrame(rows).to_csv(out/"pooled_visits.csv",index=False)

def frozen(name):
    path="b780370:results/all_species_forest/wet_network/"+name
    data=subprocess.check_output(["git","show",path],cwd=repo)
    return pd.read_csv(io.BytesIO(data),compression="gzip" if name.endswith(".gz") else None)

plants=frozen("species.csv")
flowers=frozen("flower_sources.csv.gz")
diet=frozen("diet_sources.csv.gz")
d=diet[diet.guild.eq("spiders") & diet.prey_stage.str.lower().isin(["adult","adults","adult(s)","imago"])
       & diet.dataset.notna() & ~diet.dataset.str.lower().str.contains("trophich",na=False)]
links=flowers[flowers.wet_site].merge(d,on="visitor").merge(plants[["accepted","candidate"]],left_on="plant",right_on="accepted")
counts=[]
for label,g in [("all_spiders",links),("Misumena_vatia",links[links.predator.eq("Misumena vatia")]),
                ("Tetragnatha_genus",links[links.predator.str.startswith("Tetragnatha ")]),
                ("aquatic_hoverfly_Misumena",links[links.aquatic_hoverfly & links.predator.eq("Misumena vatia")])]:
    for status in [True,False]:
        sub=g[g.candidate.eq(status)]
        counts.append(dict(group=label,candidate=status,plants=sub.plant.nunique(),visitors=sub.visitor.nunique(),
                           wet_flower_plants=int(plants[plants.candidate.eq(status)].wet_flower.sum())))
pd.DataFrame(counts).to_csv(out/"wet_links.csv",index=False)
