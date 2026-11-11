"""Find identical WCVP plant partners for named adult prey and their island-documented predators."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

from alan_dual80_prepare import resolve
from alan_gecko_primary_fiveyear import ROOT, ID
from alan_gecko_moth_species_edges import SOURCE_ROWS, SOURCE_SHA
from alan_globi_v011_same_fiveyear import taxonomy
sys.path.insert(0,str(ROOT/"scripts/acquisition"))
from alan_dual80_occurrence import Geography, SCALES, chunks, query_cells


def triples(inputs, private, edges):
    receipt=json.loads((edges/"receipt.json").read_text())
    assert receipt["scanned"]==SOURCE_ROWS and receipt["source_sha256"]==SOURCE_SHA
    for name,key in (("predator_flowers.csv","predator_flower_sha256"),("species_prey.csv","prey_table_sha256")):
        assert hashlib.sha256((edges/name).read_bytes()).hexdigest()==receipt[key]
    feeding=pd.read_csv(edges/"species_prey.csv",keep_default_na=False)
    flowers=pd.read_csv(edges/"predator_flowers.csv")
    cohort=set(pd.read_csv(ROOT/"results/dual_role_80/cohort.csv").animal_species)
    assert len(cohort)==80 and set(feeding.predator)==cohort
    direct=flowers[flowers.role.eq("reported_pollinates")][["predator","plant"]].drop_duplicates()
    assert set(direct.predator)==cohort
    adult=feeding.loc[feeding.source_stage.eq("adult") & feeding.moth.ne(""),["predator","moth"]].drop_duplicates()
    legacy=pd.read_parquet(private/"visits.parquet")
    visits=adult.merge(legacy.rename(columns={"insect_name":"moth"}),on="moth")
    assert hashlib.sha256((inputs/"wcvp.zip").read_bytes()).hexdigest()=="d32ea2b3a85e489b14e83bcc9eae7274532e1d113753f7be290d4b2dfde573fa"
    w=taxonomy(inputs,set(visits.plant_name),set(direct.plant))
    w.to_csv(private/"plant_names.csv")
    ids=w.loc[w.resolution.eq("resolved"),"accepted_id"]
    visits[ID]=visits.plant_name.map(ids)
    direct[ID]=direct.plant.map(ids)
    common=visits.dropna(subset=[ID]).merge(direct.dropna(subset=[ID]),on=["predator",ID]).drop_duplicates()
    common["accepted_name"]=common[ID].map(w[w.resolution.eq("resolved")].drop_duplicates("accepted_id").set_index("accepted_id").accepted_name)
    plants=pd.read_parquet(private/"plants.parquet")
    plants["in_primary_numerator"] = plants.wcvp_accepted_name.isin(legacy.plant_name)
    possible=common[["predator","moth",ID,"accepted_name"]].drop_duplicates().merge(plants[["key",ID,"in_primary_numerator"]],on=ID)
    return adult,visits,common,possible,receipt


def occurrence(frame,role,scale,taxa,private,geo,refs):
    dims=list(SCALES[scale])
    keys=[role,"year","land_id",*dims]
    if frame.empty:return frame.copy()
    jobs=[]
    for name,g in frame.groupby(role):
        t=taxa[name]
        if "key" not in t:raise ValueError(f"Unresolved species {name}: {t}")
        for y,geom in chunks(g,scale):
            jobs.append((y,geom,t))
    def search(job):
        y,geom,t=job
        hits,r=query_cells(private/"pages",geo,y,scale,geom,{t["key"]:t},False,adult=(role=="moth"),cell_presence=True)
        return list(hits),r
    matches=[]
    with ThreadPoolExecutor(max_workers=6) as pool:
        for data,tokens in pool.map(search,jobs):
            matches.extend(data);refs.extend(tokens)
    events=pd.DataFrame(matches,columns=keys).drop_duplicates()
    if events.empty:return frame.iloc[:0].copy()
    return frame.merge(events,on=keys,validate="many_to_one")


def counts(d):
    return {side:int(d.loc[d.land_class.eq(side),"key"].nunique()) for side in ("island","mainland")}


def run(inputs,private,edges,out):
    screen=json.loads((private/"island_screen.json").read_text())
    fine,wide=(set(screen["positives"][scale]) for scale in ("0.05deg","1deg"))
    assert fine<=wide
    adult,visits,common,possible,source=triples(inputs,private,edges)
    sites=pd.read_parquet(private/"sites.parquet")
    digest=hashlib.sha256(("\n".join(sorted(sites.target_id))+"\n").encode()).hexdigest()
    assert digest=="37f3c2c519871b6905ffa81c5a14536b248611861ad143b1d4d77d95b7c6b00d"
    base=possible[possible.predator.isin(wide)].merge(sites,on="key")
    print("Before geographic confirmation:",len(base),"rows,",base[["predator","moth","accepted_name"]].drop_duplicates().shape[0],"triples",flush=True)
    taxa={t["name"]:t for t in json.loads((private/"taxa.json").read_text())}
    for name in sorted(set(base.moth)-taxa.keys()):
        taxa[name]=resolve(name)
    geo=Geography(inputs/"coast")
    refs=[]
    stages=[]
    eligible=base
    groups=[]
    for scale in ("1deg","0.05deg"):
        if scale=="0.05deg":eligible=eligible[eligible.predator.isin(fine)]
        moths=occurrence(eligible,"moth",scale,taxa,private,geo,refs)
        detected=occurrence(moths,"predator",scale,taxa,private,geo,refs)
        print(scale,"sites",counts(detected),flush=True)
        for f,level in ((eligible,"plant"),(moths,"adult_moth_and_plant"),(detected,"database_roles_including_self_links")):
            stages.append(dict(scale=scale,level=level,site_years=counts(f),site_years_fine_cohort=counts(f[f.predator.isin(fine)])))
        groups.append(detected.assign(scale=scale))
        eligible=detected
    m=pd.concat(groups,ignore_index=True)
    m.to_parquet(private/"matched.parquet",index=False)
    review=pd.read_csv(ROOT/"results/dual_role_80/source_review.csv")
    reviewed=common[common.predator.isin(wide)].merge(review,on=["predator","moth"],how="left",validate="many_to_one")
    assert reviewed.assessment.notna().all()
    m=m.merge(review,on=["predator","moth"],validate="many_to_one")
    valid=m[m.predator.ne(m.moth) & m.assessment.eq("predation_recorded")]
    checked_totals={}
    for scale in SCALES:
        g=valid[valid.scale.eq(scale)]
        checked_totals[scale]=dict(fixed_57_species=counts(g[g.predator.isin(fine)]),
                                   scale_specific_island_cohort=counts(g),
                                   years={side:sorted(g.loc[g.land_class.eq(side),"year"].unique().tolist()) for side in ("island","mainland")})
    audit=[]
    frame=pd.read_csv(ROOT/"results/dual_role_80/cohort.csv")
    for name in frame.animal_species:
        subset=valid[valid.predator.eq(name)]
        audit.append(dict(predator=name,gbif_species_key=taxa[name].get("key"),gbif_taxonomic_status=taxa[name]["status"],
            island_005=name in fine,island_1deg=name in wide,
            named_adult_prey=adult[adult.predator.eq(name)].moth.nunique(),
            adult_prey_with_visit=visits[visits.predator.eq(name)].moth.nunique(),
            shared_plants=common[common.predator.eq(name)][ID].nunique(),
            shared_plants_on_frozen_sites=possible[possible.predator.eq(name)][ID].nunique(),
            fine_sites=counts(subset[subset.scale.eq("0.05deg")]) if name in fine else None,
            coarse_sites=counts(subset[subset.scale.eq("1deg")]) if name in wide else None))
    detail=[]
    for (name,moth,plant),grp in possible[possible.predator.isin(wide)].groupby(["predator","moth","accepted_name"]):
        sub=m[(m.predator==name)&(m.moth==moth)&(m.accepted_name==plant)]
        assessment=review[(review.predator==name)&(review.moth==moth)].assessment.item()
        detail.append(dict(predator=name,moth=moth,plant=plant,source_assessment=assessment,fine=counts(sub[sub.scale=="0.05deg"]) if name in fine else None,
                           coarse=counts(sub[sub.scale=="1deg"])))
    reviewed.to_csv(out/"island_cohort_reviewed_links.csv",index=False)
    common.to_csv(out/"shared_plant_links.csv",index=False)
    w=pd.read_csv(private/"plant_names.csv",keep_default_na=False)
    resolution={role:w.loc[w[col],"resolution"].value_counts().to_dict() for role,col in (("moth_visits","old"),("predator_pollinates","new"))}
    (out/"summary.json").write_text(json.dumps(dict(frozen_site_sha256=digest,source=source,cohort_sites=counts(sites),
        source_checked_site_years=checked_totals,wcvp_resolution=resolution, island_species_005=len(fine),island_species_1deg=len(wide),stages=stages,cohort=audit,chains=detail,
        moth_taxa={name:taxa[name].get("key") for name in sorted(set(base.moth))},
        gbif_query_pages=dict(screen=len(screen["pages"]),chains=len(refs)),
        unfiltered_database_matched_years={scale:{side:sorted(m.loc[m.scale.eq(scale)&m.land_class.eq(side),"year"].unique().astype(int).tolist())
                     for side in ("island","mainland")} for scale in SCALES}),indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    (private/"chain_pages.json").write_text(json.dumps(refs))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs","private","edges","out"):p.add_argument(arg,type=Path)
    run(**vars(p.parse_args()))
