"""Match interaction records to the species roster of long-running wet-plant surveys."""
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from alan_wet_water import stations

EARTH_KM = 6371.009


def spatial_matches(events, sites):
    chunks = []
    for year, g in events.groupby("y"):
        ref = sites[sites.y.eq(year)].reset_index(drop=True)
        if ref.empty:
            continue
        tree = BallTree(np.radians(ref[["latitude","longitude"]]),metric="haversine")
        idx, dist = tree.query_radius(np.radians(g[["lat","lon"]]),10/EARTH_KM,return_distance=True)
        if not sum(map(len,idx)):
            continue
        chunk = g.iloc[np.repeat(np.arange(len(g)),[len(i) for i in idx])].copy()
        chunk["unit_id"] = ref.unit_id.to_numpy()[np.concatenate(idx)]
        chunk["km"] = np.concatenate(dist)*EARTH_KM
        chunks.append(chunk)
    if chunks:
        return pd.concat(chunks,ignore_index=True)
    return events.head(0).assign(unit_id=pd.Series(dtype=int),km=pd.Series(dtype=float))


def run(db, d, plants, work, out):
    out.mkdir(exist_ok=True)
    net = pd.read_csv(out.parent/"wet_network/flower_links.csv.gz")[["plant","visitor"]].drop_duplicates()
    sites = pd.read_csv(out.parent/"wet_water/coordinate_trends.csv.gz")[["unit_id","island","light","state"]]
    assert sites.unit_id.nunique()==len(sites)==1707
    surveydb = duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    stations(surveydb,out.parent.parent)
    water = pd.read_csv(out.parent/"wet_water/satellite_survey_years.csv.gz").query("status=='ok'")
    surveydb.execute("""create temp table validpanel as select s.* from surveys s join both_sites using(site)
      semi join water using(unit_id,y) qualify count(*) over(partition by site)>=5 and
      max(y) over(partition by site)-min(y) over(partition by site)>=20""")
    roster = surveydb.sql("""with seen as (select distinct site,sp,candidate from year_hits join wet using(sp)
      semi join validpanel using(site,y)), both_groups as (select site from seen group by site having
      bool_or(candidate) and bool_or(not candidate)) select distinct unit_id,y,sp,candidate
      from validpanel join seen using(site) semi join both_groups using(site) semi join sites using(unit_id)""").df()
    roster = roster.merge(plants[["accepted"]].rename(columns={"accepted":"plant"}),on="sp",validate="m:1")
    xy = surveydb.sql("select distinct unit_id,latitude,longitude from geo").df()
    coords = roster[["unit_id","y"]].drop_duplicates().merge(xy,validate="m:1")
    assert coords.unit_id.nunique()==1707
    # Retain visits with known year but unavailable day as eligible annual observations.
    dopi_year = d.loc[d.Year.notna() & d.event_date.isna(),["plant","visitor","Latitude","Longitude","Year","archive","wet_site"]].copy()
    dopi_year["event_date"] = pd.to_numeric(dopi_year.Year,errors="coerce").astype("Int64").astype("string")
    supplement = dopi_year.drop(columns="Year")
    visits = db.sql("""select plant,partner as visitor,try_cast(lat as double) lat,
      try_cast(lon as double) lon,event_date,dataset,wet_site from dated_visits
      union all select plant,visitor,Latitude,Longitude,event_date,archive,wet_site from supplement""").df()
    visits = visits.merge(net,on=["plant","visitor"],validate="m:1")
    diets = db.sql("""select distinct visitor,predator,guild,prey_stage,try_cast(lat as double) lat,
      try_cast(lon as double) lon,event_date,dataset from pred""").df()
    def dated(x):
        # Single ISO year/month/day; exclude ranges and coordinates outside WGS84.
        iso=x.event_date.astype("string").str.extract(r"^(\d{4})(?:-\d{2}(?:-\d{2})?)?(?:T[^/]*)?$")[0]
        x=x.assign(y=pd.to_numeric(iso,errors="coerce"))
        return x[x.y.isin(coords.y) & x.lat.abs().le(90) & x.lon.abs().le(180) & (x.lat.ne(0)|x.lon.ne(0))].drop_duplicates()
    visits, diets = dated(visits),dated(diets)
    flower = spatial_matches(visits,coords).merge(roster,on=["plant","unit_id","y"],validate="m:1")
    pred = spatial_matches(diets,coords)
    flower.to_csv(out/"nearby_visits.csv.gz",index=False)
    pred.to_csv(out/"nearby_diet.csv.gz",index=False)
    by_species = []
    for radius in [1,10]:
        for plant,g in flower[flower.candidate & flower.km.le(radius)].groupby("plant"):
            distinct = g[["unit_id","y"]].drop_duplicates().merge(sites,validate="m:1")
            by_species.append(dict(radius_km=radius,plant=plant,site_years=len(distinct),
                sites=distinct.unit_id.nunique(),visit_years=g.y.nunique(),
                sites_with_repeated_years=int(distinct.groupby("unit_id").y.nunique().gt(1).sum()),
                visitors=g.visitor.nunique(),maintained_water_site_years=int(distinct.state.eq("maintained").sum()),
                islands=int(distinct.loc[distinct.island.eq(1),"unit_id"].nunique())))
    pd.DataFrame(by_species).to_csv(out/"candidate_visit_support.csv",index=False)
    rows=[]
    joined=[]
    for radius in [1,10]:
        f=flower[flower.km.le(radius)]
        p=pred[pred.km.le(radius)]
        for scope,g in [("any_habitat",f),("confirmed_wet",f[f.wet_site.eq(True)])]:
            for candidate in [True,False]:
                subset=g[g.candidate.eq(candidate)]
                chains=subset.merge(p,on=["unit_id","y","visitor"],suffixes=("_visit","_diet"))
                if len(chains):
                    joined.append(chains.assign(scope=scope,radius_km=radius))
                rows.append(dict(radius_km=radius,scope=scope,candidate=candidate,
                    visit_species=subset.plant.nunique(),visit_site_years=len(subset[["unit_id","y"]].drop_duplicates()),
                    shared_site_years_any_prey=len(subset[["unit_id","y"]].drop_duplicates().merge(p[["unit_id","y"]].drop_duplicates())),
                    chain_species=chains.plant.nunique(),chain_site_years=len(chains[["unit_id","y"]].drop_duplicates())))
    pd.DataFrame(rows).to_csv(out/"coverage.csv",index=False)
    if joined:
        pd.concat(joined).to_csv(out/"nearby_chains.csv.gz",index=False)
    pd.DataFrame([dict(sites=len(sites),site_years=len(coords),start=int(coords.y.min()),end=int(coords.y.max()),
       named_visits_dated_located=len(visits),wet_visits_dated_located=int(visits.wet_site.eq(True).sum()),
       diet_dated_located=len(diets),diet_site_years_1km=len(pred[pred.km.le(1)][["unit_id","y"]].drop_duplicates()),
       diet_site_years_10km=len(pred[["unit_id","y"]].drop_duplicates()))]).to_csv(out/"survey_support.csv",index=False)
