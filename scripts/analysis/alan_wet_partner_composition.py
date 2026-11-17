"""Source-matched recorded interaction composition in the frozen wet cohort."""
import io
import subprocess
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from statsmodels.stats.multitest import multipletests

root = Path(__file__).resolve().parents[2]
base = root / "results/all_species_forest"
out = base / "wet_partner_composition"
out.mkdir(exist_ok=True)

def frozen(name):
    raw = subprocess.check_output(["git", "show", "b780370:results/all_species_forest/wet_network/" + name], cwd=root)
    return pd.read_csv(io.BytesIO(raw), compression="gzip" if name.endswith(".gz") else None)

def records(filename):
    df = frozen(filename).dropna(subset=["dataset"])
    return df[~df.dataset.str.contains("trophich", case=False)]

w = frozen("species.csv")[["sp","species","accepted","indicator","candidate","site_region"]]
assert len(w)==363 and w.candidate.sum()==58 and w.sp.is_unique
h = pd.read_csv(base/"seven_moth/habitat_optima.csv.gz").groupby("sp").optimum.agg(set)
flowers = records("globi_flower_sources.csv.gz").rename(columns={"dataset":"flower_source"}).merge(
    w[~w.accepted.duplicated(keep=False)][["sp","accepted"]],left_on="plant",right_on="accepted",validate="m:1")[["sp","flower_source","visitor"]].drop_duplicates()
diet = records("diet_sources.csv.gz").rename(columns={"dataset":"diet_source"})
adult = diet.prey_stage.str.lower().isin(["adult","adults","adult(s)","imago"])
pair_results = []
scores = []
indicators = []
source_rows = []
inventories = []
for layer, selected in [("visitors",None),("predators",diet),("adult_prey_predators",diet[adult])]:
    if selected is None:
        links = flowers.rename(columns={"visitor":"partner"}).assign(diet_source="")
    else:
        feeding = selected[["visitor","diet_source","predator"]].drop_duplicates()
        links = flowers.merge(feeding,on="visitor").drop(columns="visitor").rename(columns={"predator":"partner"}).drop_duplicates()
    links = links[["sp","flower_source","diet_source","partner"]].drop_duplicates()
    links["layer"] = layer
    inventories.append(links)
    counts = links.groupby("sp").partner.nunique()
    w[layer+"_documented"] = w.sp.map(counts).fillna(0).astype(int)
    supported = w[w.sp.isin(counts.index)].reset_index(drop=True)
    ix = dict(zip(supported.sp,supported.index))
    label = supported.candidate.to_numpy()
    rng = np.random.default_rng(9112)
    # Shuffle plant labels once per permutation, identically across all source panels.
    richness = supported.sp.map(counts)
    strata = supported.groupby(["site_region",pd.cut(richness,[0,4,19,np.inf],labels=False)]).indices
    labels = np.tile(label,(2000,1))
    for j in range(1,len(labels)):
        for group in strata.values():
            labels[j,group] = rng.permutation(label[group])
    pairs = {}
    for (f,d),panel in links.groupby(["flower_source","diet_source"]):
        matrix = pd.crosstab(panel.sp,panel.partner).gt(0)
        ids = np.array([ix[s] for s in matrix.index])
        ncase = label[ids].sum()
        keep = ncase>0 and ncase<len(ids)
        source_rows.append((layer,f,d,len(ids),ncase,len(ids)-ncase,len(matrix.columns),keep))
        if not keep: continue
        for (a,b),dist in zip(combinations(ids,2),pdist(matrix.to_numpy(),metric="jaccard")):
            pairs.setdefault(tuple(sorted((a,b))),[]).append(1-dist)
    # Each eligible plant pair has equal weight regardless of number of shared sources.
    arr = np.array(list(pairs))
    sim = np.array([np.mean(x) for x in pairs.values()])
    habitats = supported.sp.map(h)
    different = []
    for (a,b),v in pairs.items():
        u,z = habitats[a],habitats[b]
        distinct = isinstance(u,set) and isinstance(z,set) and u.isdisjoint(z)
        different.append(distinct)
        pair_results.append((layer,supported.sp[a],supported.sp[b],len(v),np.mean(v),distinct))
    for subset,mask in [('all',np.ones(len(sim),bool)),('nonoverlapping_optima',np.array(different))]:
        a,b=arr[mask].T
        sims=sim[mask]
        stat=[]
        for l in labels:
            cc=l[a]&l[b]
            cr=l[a]^l[b]
            stat.append(sims[cc].mean()-sims[cr].mean() if cc.any() and cr.any() else np.nan)
        cc=label[a]&label[b];cr=label[a]^label[b];rr=~(cc|cr)
        perm=np.array(stat[1:]);perm=perm[np.isfinite(perm)]
        scores.append((layer,subset,len(set(a[cc])|set(b[cc])),len(set(a[rr])|set(b[rr])),cc.sum(),cr.sum(),rr.sum(),sims[cc].mean(),sims[cr].mean(),sims[rr].mean(),stat[0],np.mean(perm),(1+sum(perm>=stat[0]-1e-12))/(1+len(perm))))
    # Describe the partners shared by candidate pairs with disjoint habitat optima.
    # Intersections are computed within source panels before deduplication across sources.
    by_partner = {}
    for _,panel in links.groupby(["flower_source","diet_source"]):
        if panel.sp.map(ix).map(lambda i:label[i]).nunique()<2:continue
        for partner,g in panel.groupby("partner"):
            ids=sorted({ix[x] for x in g.sp})
            for a,b in combinations(ids,2):
                if isinstance(habitats[a],set) and isinstance(habitats[b],set) and habitats[a].isdisjoint(habitats[b]):
                    by_partner.setdefault(partner,set()).add((a,b))
    for partner,pairset in by_partner.items():
        ccset = [(a,b) for a,b in pairset if label[a] and label[b]]
        mixed = sum(label[a] != label[b] for a,b in pairset)
        indicators.append((layer,partner,len(ccset),len({k for pair in ccset for k in pair}),mixed,len(pairset)-len(ccset)-mixed))

w["shared_accepted_name_excluded"] = w.accepted.duplicated(keep=False)
w.to_csv(out/"species.csv",index=False)
pd.concat(inventories).to_csv(out/"links.csv.gz",index=False)
pd.DataFrame(source_rows,columns="layer flower_source diet_source plants candidates references partners eligible".split()).to_csv(out/"sources.csv",index=False)
pd.DataFrame(pair_results,columns="layer first second shared_panels similarity disjoint_known_optima".split()).to_csv(out/"pairs.csv.gz",index=False)
tests=pd.DataFrame(scores,columns="layer subset candidates_in_cc references_in_rr cc_pairs mixed_pairs rr_pairs cc_similarity mixed_similarity rr_similarity difference perm_mean p".split())
tests["q"]=multipletests(tests.p,method="fdr_bh")[1]
tests.to_csv(out/"tests.csv",index=False)
pd.DataFrame(indicators,columns="layer partner candidate_pairs candidate_species mixed_pairs reference_pairs".split()).sort_values(["layer","candidate_pairs"],ascending=[True,False]).to_csv(out/"shared_partners.csv",index=False)
print(tests.to_string(index=False))
