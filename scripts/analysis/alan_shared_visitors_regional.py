"""Join effort-documented flower surveys to the frozen five-year regional comparison."""
import argparse
import geopandas as gpd
import json
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.stats import t as student_t

from alan_gecko_primary_fiveyear import primary, ID, assign_points, prepare_land
from alan_globi_v011_same_fiveyear import taxonomy

YEARS={2019,2021,2023,2024,2025}
HDC=["Hymenoptera","Diptera","Coleoptera"]
TRAITS=["logpop","logphotos","logspecies","bio5","lograin","elevation_km"]


def grid(x,y):
    return (np.floor((x+180)/5).astype(int).astype(str)+"_"+np.floor((y+90)/5).astype(int).astype(str))


def uk_coordinates(ref,proj):
    letters="ABCDEFGHJKLMNOPQRSTUVWXYZ"
    m=re.fullmatch(r"([A-Z]{1,2})(\d{2})(\d{2})",str(ref))
    if not m:return np.nan,np.nan
    code,e,n=m.groups()
    if proj=="OSGB" and len(code)==2:
        i,j=map(letters.index,code)
        east=((i-2)%5)*5 + j%5
        north=19-(i//5)*5-j//5
        crs=27700
    elif proj=="OSIE" and len(code)==1:
        j=letters.index(code); east=j%5; north=4-j//5;crs=29903
    else:return np.nan,np.nan
    return TRANS[crs].transform(east*100000+int(e)*1000+500,north*100000+int(n)*1000+500)


TRANS={c:Transformer.from_crs(c,4326,always_xy=True) for c in [27700,29903]}


def uk(p):
    with zipfile.ZipFile(p/'uk.zip') as z:
        d=pd.concat([pd.read_csv(z.open(n),encoding="cp1252") for n in z.namelist() if n.endswith(".csv")],ignore_index=True)
    d=d[d.year.isin(YEARS)].copy()
    assert d.sample_id.is_unique
    target=d.target_flower_corrected.fillna("")
    label=target.str.split(" - ",n=1).str[-1]
    other=d.target_other_name_corrected.fillna(d.target_other_name).fillna("")
    d["original"]=label.where(label.str.fullmatch(r"[A-Z][a-z]+ [a-z-]+"),other)
    # UK PoMS records bramble and dandelion as aggregates.
    d.loc[target.str.startswith(("Bramble","Dandelion")),"original"]=""
    coords=d[["x1km_square","sample_projection"]].drop_duplicates().copy()
    coords[["x","y"]]=[uk_coordinates(a,b) for a,b in coords.itertuples(index=False,name=None)]
    d=d.merge(coords,on=["x1km_square","sample_projection"],validate="many_to_one")
    d["Hymenoptera"]=d[["bumblebees","honeybees","solitary_bees","wasps"]].sum(axis=1)
    d["Diptera"]=d[["hoverflies","other_flies"]].sum(axis=1);d["Coleoptera"]=d.beetles
    d["total"]=d.all_insects_total
    assert (d[HDC].sum(axis=1)<=d.total).all()
    d["source"]="UK_FIT";d["event"]=d.sample_id.astype(str);d["comparable"]=True
    return d


def european(p):
    d=pd.read_csv(p/"Interaction_data.csv.gz",encoding="latin1",low_memory=False)
    d["year"]=pd.to_numeric(d.Date.str[:4],errors="coerce")
    # Confirm record dates against the original working metadata, including direct survey rounds.
    permitted={"5_Marini":[2019],"10_Vanbergen":[2021],"11_Clough":[2021],"13_Karise":[2019,2021],
               "27_Castro":[2021],"39_Schweiger":[2021],"42_Knight":[2019],"43_Knight":[2019],
               "49_Hervias-Parejo":[2019],"50_Hervias-Parejo":[2021]}
    d=d[pd.concat([d.Study_id.eq(k)&d.year.isin(v) for k,v in permitted.items()],axis=1).any(axis=1)].copy()
    d=d[d.Interaction.gt(0) & d.Plant_unsure_id.astype(str).str.casefold().eq("no")].copy()
    d=d[d.Plant_rank.str.lower().eq("species")].copy()
    d=d.rename(columns={"Latitude":"y","Longitude":"x","Plant_accepted_name":"original","Interaction":"total"})
    for h in HDC:d[h]=d.total.where(d.Pollinator_order.eq(h),0)
    d["source"]="EuPPollNet";d["event"]=d.Study_id+"|"+d.Network_id+"|"+d.Date
    # Species/group-restricted studies are valid evidence of presence but not total visitor diversity.
    d["comparable"]=d.Study_id.isin(["5_Marini","10_Vanbergen","11_Clough","27_Castro","39_Schweiger","49_Hervias-Parejo","50_Hervias-Parejo"])
    return d


def pollagri(p):
    keep=["ID","CROP_LB_NOM","CROP_LATITUDE","CROP_LONGITUDE","STUDY_YEAR","STUDY_SAMPLING_DURATION","POLL_ORDER","SOURCE","STUDY_METHOD"]
    d=pd.read_csv(p/"fr_pollAgri.csv",sep="\t",usecols=keep,low_memory=False)
    for k in ["STUDY_SAMPLING_DURATION","CROP_LATITUDE","CROP_LONGITUDE"]:
        d[k]=pd.to_numeric(d[k].astype(str).str.replace(",","."),errors="coerce")
    d=d[d.STUDY_YEAR.isin(YEARS)&d.STUDY_SAMPLING_DURATION.gt(0)&~d.SOURCE.str.contains("Allen-Perkins",na=False)].copy()
    d=d.rename(columns={"CROP_LB_NOM":"original","STUDY_YEAR":"year","CROP_LATITUDE":"y","CROP_LONGITUDE":"x"})
    for h in HDC:d[h]=d.POLL_ORDER.eq(h).astype(int)
    d["source"]="PollAgri";d["event"]=None;d["total"]=1
    # A taxon-interaction entry does not carry a standardized number of individual visits.
    d["comparable"]=False
    return d


def crops(p):
    field=pd.read_csv(p/"CropPol_field_level_data.csv",encoding="latin1",low_memory=False)
    data=pd.read_csv(p/"CropPol_sampling_data.csv",encoding="latin1",low_memory=False)
    field["year"]=pd.to_numeric(field.sampling_year,errors="coerce")
    field=field[field.year.isin(YEARS)].copy()
    d=data.merge(field[["study_id","site_id","crop","year","latitude","longitude","use_visits_or_abundance"]],on=["study_id","site_id"])
    amount=d[["total_sampled_area","total_sampled_time","total_sampled_flowers"]].apply(pd.to_numeric,errors="coerce")
    d=d[amount.gt(0).any(axis=1)&pd.to_numeric(d.abundance,errors="coerce").gt(0)&d.use_visits_or_abundance.astype(str).str.contains("visit",case=False)].copy()
    d=d.rename(columns={"crop":"original","latitude":"y","longitude":"x"})
    d["total"]=pd.to_numeric(d.abundance);d["source"]="CropPol";d["event"]=d.study_id+"|"+d.site_id.astype(str)
    d["comparable"]=False
    for h in HDC:d[h]=0
    return d


def oregon(p):
    d=pd.read_csv(p/"sa026.csv",low_memory=False)
    d=d[d.YEAR.isin(YEARS)&d.PLTSP_NAME.notna()&d.VISSP_NAME.notna() & pd.to_numeric(d.NO_INT,errors="coerce").gt(0)].copy()
    # All study meadows are in the same 5-degree region (Andrews Forest, Oregon).
    d=d.rename(columns={"PLTSP_NAME":"original","YEAR":"year","NO_INT":"total"})
    d["x"],d["y"]=-122.2,44.2
    d["source"]="Andrews";d["event"]=d.PLOT_ID+"|"+d.SAMPLEDATE+"|"+d.WATCH.astype(str)
    d["comparable"]=False
    for h in HDC:d[h]=0
    return d


def metric(n,parts):
    if n<2:return np.nan
    return sum(1-(n-a)*(n-a-1)/(n*(n-1)) for a in parts)


def fit(frame):
    g=frame.copy()
    g["light_x_index"]=g.loglight*g.index_value
    g["island_x_index"]=g.island*g.index_value
    g["light_island_x_index"]=g.light_island*g.index_value
    cols=["island","loglight","light_island","light_x_index","island_x_index","light_island_x_index",*TRAITS]
    sizes=g.groupby(["region_year","land_class"]).size()
    w=sizes.groupby(level=0).min().reindex(g.region_year).to_numpy()/sizes.reindex(pd.MultiIndex.from_frame(g[["region_year","land_class"]])).to_numpy()
    z=g[["response",*cols]].astype(float)
    w=pd.Series(w,index=g.index)
    mu=z.mul(w,axis=0).groupby(g.region_year).transform("sum").div(w.groupby(g.region_year).transform("sum"),axis=0)
    z=(z-mu).mul(np.sqrt(w),axis=0)
    X=z[cols].to_numpy();y=z.response.to_numpy()
    rank=np.linalg.matrix_rank(X)
    info=dict(n=len(g),regions=g.region5.nunique(),region_years=g.region_year.nunique(),islands=int(g.island.sum()),mainland=int((1-g.island).sum()),rank=int(rank),parameters=len(cols))
    if rank<len(cols):return info
    beta=np.linalg.lstsq(X,y,rcond=None)[0];delta=beta[5]
    low,high=g[["region_year","index_value"]].drop_duplicates().index_value.quantile([.25,.75])
    info.update(interaction=delta,q25=low,q75=high,gap_q25=beta[2]+delta*low,gap_q75=beta[2]+delta*high)
    if info["regions"]>=10 and len(g)>=50:
        u=X*(y-X@beta)[:,None]
        sums=pd.DataFrame(u).groupby(g.region5.to_numpy()).sum().to_numpy()
        bread=np.linalg.inv(X.T@X)
        G=info["regions"];df=len(g)-rank-info["region_years"]
        V=bread@(sums.T@sums)@bread*(G/(G-1))*((len(g)-1)/df)
        se=np.sqrt(V[5,5]);crit=student_t.ppf(.975,G-1)
        info["ci95"]=[delta-crit*se,delta+crit*se]
    return info


def run(inputs,signed,p,out):
    out.mkdir(parents=True,exist_ok=True)
    s,digest,plants,legacy=primary(inputs,signed)
    plants=plants.merge(s[["key","land_class","year","region5"]],on="key",validate="many_to_one")
    yes=plants[plants.all_moth]
    shared=yes.groupby("wcvp_accepted_name").land_class.nunique().loc[lambda v:v==2].index
    assert len(shared)==544
    yes=yes[yes.wcvp_accepted_name.isin(shared)]
    s["response"]=100*s.key.map(yes.groupby("key")[ID].nunique()).fillna(0)/s.n_species
    f=pd.concat([d[["source","event","year","x","y","original","total",*HDC,"comparable"]] for d in (uk(p),european(p),pollagri(p),crops(p),oregon(p))],ignore_index=True)
    f=f.dropna(subset=["x","y","original"]);f["region5"]=grid(f.x,f.y)
    names={n.casefold():n for n in shared}
    # Resolve species-level Latin labels only; no genus/group/common-name guesses.
    f["latin"]=f.original.astype(str).str.strip().str.lower().str.fullmatch(r"[a-z]+ [a-z-]+")
    raw=f.loc[f.latin,"original"].str.strip().str.lower().unique()
    canonical={n:n.capitalize() for n in raw}
    w=taxonomy(inputs,set(canonical.values()),set())
    resolved={k:v.accepted_name for k,v in w.iterrows() if v.resolution=="resolved"}
    labels={k:names.get(k,resolved.get(v)) for k,v in canonical.items()}
    f["plant"]=f.original.astype(str).str.strip().str.lower().map(labels)
    f=f[f.plant.isin(shared)].copy()
    keys=["year","region5","plant"]
    local=yes[["year","region5","wcvp_accepted_name","land_class"]].drop_duplicates().rename(columns={"wcvp_accepted_name":"plant"})
    g=f.merge(local[keys].drop_duplicates(),on=keys)
    def extent(d):return dict(species=d.plant.nunique(),regions=d.region5.nunique(),region_years=d[["region5","year"]].drop_duplicates().shape[0],species_region_years=d[keys].drop_duplicates().shape[0])
    # All timed counts contribute effort, including a zero-visitor session.
    coverage=extent(g[g.total.gt(0)]);coverage["by_source"]={k:extent(v[v.total.gt(0)]) for k,v in g.groupby("source")}
    allcov={k:extent(v[v.total.gt(0)]) for k,v in f.groupby("source")}
    positive=g[g.total.gt(0)][keys+["source"]].drop_duplicates()
    positive.to_csv(out/"regional_plant_support.csv",index=False)
    species=pd.DataFrame({"plant":shared});ns=positive.groupby("plant").region5.nunique()
    species["matched_regions"]=species.plant.map(ns).fillna(0).astype(int);species.to_csv(out/"shared_544.csv",index=False)
    comp=g[g.comparable];by=comp.groupby(["source",*keys])[["total",*HDC]].sum().reset_index()
    by["diversity"]=[metric(n,a) for n,a in zip(by.total,by[HDC].to_numpy())]
    by=by.dropna(subset=["diversity"])
    def regional(data,col):
        d=data.groupby(keys)[col].mean().groupby(["year","region5"]).mean().rename("index_value").reset_index()
        return d
    e=regional(by,"diversity")
    fit_results={}
    def analyze(d,label):
        frame=s.merge(d,on=["region5","year"],validate="many_to_one")
        if frame.empty:return {}
        d=d.copy()
        ranges=frame.groupby(["region_year","land_class"]).loglight.agg(["min","max","nunique"]).unstack()
        # Slopes within each land class require variation, and the comparison should share light support.
        valid=ranges.index[(ranges["nunique"].min(axis=1)>=2)&(ranges["min"].max(axis=1)<ranges["max"].min(axis=1))]
        supported=frame[frame.region_year.isin(valid)]
        res={"all_matched":fit(frame),"within_land_light_overlap":fit(supported) if len(supported) else None}
        d["light_overlap"]=d.year.astype(int).astype(str).str.cat(d.region5,sep="|").isin(valid)
        d.to_csv(out/(label+"_regions.csv"),index=False)
        return res
    fit_results["order_diversity"]=analyze(e,"order_diversity")
    ukg=comp[comp.source.eq("UK_FIT")].copy()
    ukg["non_lep_visitors_per_10min"]=ukg[HDC].sum(axis=1)
    rate=regional(ukg,"non_lep_visitors_per_10min")
    fit_results["uk_visitors_per_10min"]=analyze(rate,"uk_frequency")
    s["response"]=100*s.all_moth/s.n_species
    fit_results["original_response_order_diversity"]=analyze(e,"order_diversity")
    fit_results["original_response_uk_frequency"]=analyze(rate,"uk_frequency")
    both=local.groupby(keys).land_class.nunique().loc[lambda x:x==2].reset_index()[keys]
    coords=g[["x","y"]].drop_duplicates().reset_index(drop=True).rename(columns={"x":"longitude","y":"latitude"})
    coords["unit_id"]=coords.index
    coast=inputs/"coast"
    land=assign_points(coords,prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp")),gpd.read_file(coast/"GSHHS_h_L2.shp"))
    survey=g.merge(coords.merge(land,on="unit_id").rename(columns={"longitude":"x","latitude":"y"})[["x","y","land_class"]],on=["x","y"])
    survey=survey[survey.total.gt(0)]
    survey[["year","region5","source","plant","land_class"]].drop_duplicates().to_csv(out/"survey_land_support.csv",index=False)
    two=survey[survey.land_class.isin(["island","mainland"])].groupby(["region5","year"]).land_class.nunique()
    lightok=pd.read_csv(out/"order_diversity_regions.csv").query("light_overlap")[["year","region5"]]
    summary=dict(target_sha256=digest,condition="same plant, same 5-degree region and observed year",survey_coverage_before_cohort=allcov,
                 matched_effort_including_zero=extent(g),matched_positive_coverage=coverage,local_photo_both_sides=extent(positive.merge(both,on=keys)),survey_region_years_both_sides=int(two.eq(2).sum()),comparable_coverage=extent(by),matched_light_overlap=extent(by.merge(lightok,on=["year","region5"])),
                 uk_timed_counts=len(ukg),uk_minutes=int(10*len(ukg)),support_by_year=positive.groupby("year").size().to_dict(),models=fit_results)
    (out/"summary.json").write_text(json.dumps(summary,indent=2,default=lambda x:x.item())+"\n")
    print(json.dumps(summary,default=lambda x:x.item()),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ["inputs","signed","p","out"]:p.add_argument(k,type=Path)
    run(**vars(p.parse_args()))
