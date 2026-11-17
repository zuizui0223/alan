"""Compare sites paired on light trend and habitat, independently of all outcomes."""
import networkx as nx
import numpy as np
import pandas as pd
from alan_wet_shore import repo
out=repo/"results/all_species_forest/wet_annual"
s=pd.read_csv(out/"site_trends.csv").set_index("unit_id")
edges=nx.Graph()
for i,a in s.iterrows():
    for j,b in s.loc[s.index>i].iterrows():
        cost=abs(a.light-b.light)
        if a.habitat==b.habitat and cost<=1 and np.sign(round(a.light,8))==np.sign(round(b.light,8)):
            edges.add_edge(i,j,weight=2-cost)
pairs=sorted(tuple(sorted(x)) for x in nx.max_weight_matching(edges,maxcardinality=True))
rows=[]
for i,j in pairs:
    for state in ["drift_water_freq","drift_ndvi","turnover"]:
        a,b=s.loc[i],s.loc[j]
        if b[state]>a[state]:a,b=b,a
        if not np.isfinite(a[state]-b[state]) or a[state]-b[state]<1e-9:continue
        rows.append(dict(state=state,more=int(a.name),less=int(b.name),more_light=a.light,less_light=b.light,state_gap=a[state]-b[state],candidate_change=a.candidate_rate-b.candidate_rate,reference_change=a.reference_rate-b.reference_rate,difference=a.difference-b.difference))
d=pd.DataFrame(rows)
d.to_csv(out/"light_matched_pairs.csv",index=False)
print(d.groupby("state").difference.agg(["count","median",lambda v:(v<0).sum()]))
