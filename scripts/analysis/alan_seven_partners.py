"""Pairwise similarity of documented partners; primary pollen spectrum versus pooled GloBI."""
from itertools import combinations
from pathlib import Path
import io,subprocess
import numpy as np
import pandas as pd

repo=Path(__file__).resolve().parents[2]
root=repo/"results/all_species_forest"
out=root/"seven_moth"
w=pd.read_csv(root/"wet_reproduction/species.csv")
seven=w[w.candidate & w.documented_transport].species.tolist()
s=pd.read_excel(repo.parent/"gbif_annual/walton_rsbl20190877_si_001.xlsx",sheet_name="Plant-Pollinator Table (S3)",header=3)
s.index=s.iloc[:,0].str.strip().str.rstrip("*").replace({"Iris psuedacorus":"Iris pseudacorus"})
families=s.iloc[:,5:12].loc[seven].copy().rename(columns={"Hepliaidae":"Hepialidae"})
assert np.array_equal(families.sum(axis=1),w.set_index("species").loc[seven,"walton_s3_moths"])
families.to_csv(out/"pollen_carriers.csv")
frac=families.div(families.sum(axis=1),axis=0)
g=pd.read_csv(io.BytesIO(subprocess.check_output(["git","show","b780370:results/all_species_forest/wet_network/globi_flower_sources.csv.gz"],cwd=repo)),compression="gzip")
# Count documented insect species only within namespaces shared by the two plants.
g=g[g.plant.isin(seven)].dropna(subset=["dataset","visitor"])
g=g[~g.dataset.str.contains("trophich",case=False)]
rows=[]
for a,b in combinations(seven,2):
    x=g[g.plant.eq(a)];y=g[g.plant.eq(b)]
    common=set(x.dataset)&set(y.dataset)
    overlaps=set()
    union=set()
    for dataset in common:
        u=set(x.loc[x.dataset.eq(dataset),"visitor"])
        v=set(y.loc[y.dataset.eq(dataset),"visitor"])
        overlaps|=u&v
        union|=u|v
    rows.append(dict(first=a,second=b,moth_family_bray_curtis=abs(frac.loc[a]-frac.loc[b]).sum()/2,
        shared_namespaces=len(common),shared_visitors=len(overlaps) if common else None,
        recorded_visitor_union=len(union) if common else None))
pd.DataFrame(rows).to_csv(out/"partner_similarity.csv",index=False)
# Secondary description of the wet habitat classes (nonexclusive published optima).
from alan_rf_decline_traits import compare
from statsmodels.stats.multitest import multipletests
h=pd.read_csv(out/"habitat_optima.csv.gz")
env=pd.read_csv(root/"decline_traits/moisture_evidence.csv.gz").set_index("sp").loc[w.sp]
env["seven"]=env.species.isin(seven)
has=env.index.isin(h.sp)
habrows=[]
for name,g in h.groupby("optimum"):
    v=pd.Series(env.index.isin(g.sp),index=env.index).astype(float).where(has)
    test=compare(env,v,["site_region","island_band","site_band"],np.random.default_rng(714))
    habrows.append(dict(habitat=name,candidate=int(v[env.candidate].sum()),reference=int(v[~env.candidate].sum()),
                        seven=int(v[env.seven].sum()),candidate_known=int((env.candidate & has).sum()),
                        reference_known=int((~env.candidate & has).sum()),**test))
hab=pd.DataFrame(habrows)
hab["q"]=multipletests(hab.p,method="fdr_bh")[1]
hab.sort_values(["seven","candidate"],ascending=False).to_csv(out/"habitat_prevalence.csv",index=False)
