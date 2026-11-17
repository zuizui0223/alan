"""Identical adult-visitor and predator pairs across the two wet plant groups."""
from pathlib import Path
import io
import subprocess
import pandas as pd

repo=Path(__file__).resolve().parents[2]
base=repo/"results/all_species_forest"
out=base/"wet_shore"
w=pd.read_csv(base/"wet_partner_composition/species.csv")
vis=pd.read_csv(base/"wet_partner_composition/links.csv.gz").query("layer=='visitors'").rename(columns={"partner":"visitor"})
raw=subprocess.check_output(["git","show","b780370:results/all_species_forest/wet_network/diet_sources.csv.gz"],cwd=repo)
diet=pd.read_csv(io.BytesIO(raw),compression="gzip").dropna(subset=["dataset"])
diet=diet[diet.prey_stage.str.lower().isin(["adult","adults","adult(s)","imago"]) & ~diet.dataset.str.contains("trophich",case=False)]
diet=diet.rename(columns={"dataset":"diet_source"})
chain=vis[["sp","flower_source","visitor"]].merge(diet[["visitor","predator","guild","diet_source"]].drop_duplicates(),on="visitor").merge(w[["sp","candidate"]],on="sp")
hab=pd.read_csv(base/"seven_moth/habitat_optima.csv.gz").groupby("sp").optimum.agg(set).to_dict()
rows=[]
for (v,p,g),data in chain.groupby(["visitor","predator","guild"]):
    cset=set(); rset=set(); edges=set(); panels=0; independent=0
    for (flower,d),s in data.groupby(["flower_source","diet_source"]):
        c=set(s[s.candidate].sp); r=set(s[~s.candidate].sp)
        if not (c and r):continue
        cset.update(c);rset.update(r);panels+=1
        independent+=int("inaturalist" not in (flower+d).lower())
        edges.update((i,j) for i in c for j in r if i in hab and j in hab and hab[i].isdisjoint(hab[j]))
    if panels:rows.append((v,p,g,len(cset),len(rset),panels,independent,len(edges)))
cols="visitor predator guild candidates references common_source_pairs source_pairs_without_inat cross_group_disjoint_habitat_pairs".split()
res=pd.DataFrame(rows,columns=cols).sort_values(["candidates","references","cross_group_disjoint_habitat_pairs"],ascending=False)
res.to_csv(out/"shared_adult_chains.csv",index=False)
chain.to_csv(out/"plant_chain_records.csv.gz",index=False)
print(len(res),sum(res.cross_group_disjoint_habitat_pairs>0))
print(res.head(12).to_string(index=False))
print(res.query("source_pairs_without_inat>0").head(4).to_string(index=False))
