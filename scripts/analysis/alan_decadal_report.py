from alan_decadal_events import *
from alan_review_satellite_spring import greenup_date
import json
import statsmodels.formula.api as sm
from collections import defaultdict
from scipy.stats import theilslopes

def event_panel():
    records=[]
    for kind in ["plant","insect","bird","fruit"]:
        for f in sorted(C.glob(f"npn_*_{kind}.json")):
            records.extend(dict(kind=kind,station=r["site_id"],id=r["species_id"],individual=r["individual_id"],
                date=r["observation_date"],status=r["phenophase_status"],phase=r["phenophase_id"]) for r in json.loads(f.read_bytes()))
    a=pd.DataFrame(records).drop_duplicates()
    a=a[a.status.isin([0,1])].copy()
    a["year"]=pd.to_datetime(a.date).dt.year
    a["day"]=pd.to_datetime(a.date).dt.dayofyear
    a["taxon"]=a.id.map(CAT.set_index("species_id").taxon)
    a.loc[a.kind.isin(["plant","fruit"]),"taxon"]=a.id.map(names.set_index("species_id").accepted_name)
    a=a.merge(coords.land_class.rename("land"),left_on="station",right_index=True)
    on=[]
    for (st,t,y,ind),g in a[a.kind.eq("plant")].groupby(["station","taxon","year","individual"]):
        if t=="Hamamelis virginiana":g=g[g.day>=183]
        first=g.loc[g.status.eq(1),"day"].min()
        last=g.loc[g.status.eq(0)&g.day.lt(first),"day"].max()
        if first-last<=30:on.append([st,t,y,(first+last)/2])
    onset=pd.DataFrame(on,columns=["station","plant","year","onset"]).groupby(["station","plant","year"],as_index=False).onset.median()
    links=pd.read_csv(OUT/"bird_links.csv")
    plantlinks=LEG[['plant_name','insect_name','guild']].drop_duplicates().set_axis(['plant','animal','guild'],axis=1)
    birds=links[['target_taxon_name','source_taxon_name','interaction_type']].set_axis(['plant','animal','guild'],axis=1)
    links=pd.concat([plantlinks,birds[birds.guild.isin(['eats','visitsFlowersOf','pollinates'])]]).drop_duplicates()
    dates=a[a.status.eq(1)].groupby(['kind','station','taxon','year']).day.agg(set).to_dict()
    rows=[]
    for t in targets.itertuples():
        for link in links[links.plant.eq(t.plant)].itertuples():
            animal='insect' if 'family_visit' in link.guild else 'bird'
            for y in range(2008,2025):
                plantdays=dates.get(('fruit' if link.guild=='eats' else 'plant',t.station,t.plant,y),set())
                days=dates.get((animal,t.station,link.animal,y),set())
                if days and plantdays:
                    rows.append([t.station,t.land_class,t.plant,link.animal,link.guild,y,len(plantdays),len(days),
                        len(plantdays & days),np.median(list(days))-np.median(list(plantdays))])
    overlaps=pd.DataFrame(rows,columns=['station','land','plant','animal','guild','year','plant_dates','animal_dates','same_dates','animal_minus_plant_median_observed_day'])
    overlaps.to_csv(OUT/'species_date_overlap.csv',index=False)
    return onset

def fit(d,field,unit):
    d=d.dropna(subset=[field,'light_dn']).copy()
    d=d[d.groupby(unit).year.transform('size')>=2].copy()
    d['u']=d[unit].astype(str);d['yr']=d.year.astype(str)
    d['x']=d.light_dn/10; d['ix']=d.x*d.land.eq('island')
    f=sm.ols(f'{field} ~ x + ix + C(u) + C(yr)',d,eval_env=-1).fit(cov_type='cluster',cov_kwds={'groups':d.station},use_t=True)
    out=[]
    for side,term in [('mainland','x'),('island','x+ix'),('difference','ix')]:
        t=f.t_test(term+'=0');out.append(dict(metric=field,side=side,estimate=t.effect.item(),lo=t.conf_int()[0,0],hi=t.conf_int()[0,1],n=len(d),stations=d.station.nunique(),first=d.year.min(),last=d.year.max()))
    return out

def run():
    light=pd.read_csv(OUT/'light.csv')
    env=light.copy()
    green={}
    manifest=OUT/'green_source_batches.csv'
    paths=[C/name for name in pd.read_csv(manifest).file] if manifest.exists() else sorted(C.glob('green_*.json'))
    print('Retrieved greenup batches',len(paths),'of',3*len(coords),flush=True)
    for path in paths:
        _,st,_=path.stem.split('_');obj=json.loads(path.read_bytes())
        assert abs(obj['latitude']-coords.loc[int(st),'latitude'])<1e-5 and abs(obj['longitude']-coords.loc[int(st),'longitude'])<1e-5
        years=defaultdict(dict)
        for r in obj['subset']:years[int(r['modis_date'][1:5])][r['band']]=r['data'][0]
        for y,bands in years.items():
            day=greenup_date(bands)
            if pd.notna(day) and day.year==y:green[int(st),y]=day.dayofyear
    summer=defaultdict(dict)
    for path in sorted((F/"extended11").glob("px_*.json"))+sorted((F/"long_ndvi").glob("ndvi_*.json")):
        station=int(path.name.split("_")[1]); year=int(path.name.split(".")[1][1:5])
        v=json.loads(path.read_bytes())["values"]
        raw,qa,day=v
        if qa==0 and day and 160<=day<=240 and raw is not None and -2000<=raw<=10000:
            summer[station,year][day]=raw/10000
    env["summer_ndvi"]=[np.median(list(summer[k].values())) if len(summer[k])>=2 else None for k in zip(env.station,env.year)]
    env['green_doy']=[green.get((r.station,r.year)) for r in env.itertuples()]
    def biomass(s,y):
        f=C/f'woody_{s}_{y}.json'
        if not f.exists():return None
        v=json.loads(f.read_bytes())['values'][0]
        return v if v is not None and 0<=v<2147483647 else None
    env['woody_tonnes']=[biomass(r.station,r.year) for r in env.itertuples()]
    # A common grid cell must not count as multiple independent biomass observations.
    import math
    meta=requests.get('https://planetarycomputer.microsoft.com/api/stac/v1/collections/chloris-biomass/items/chloris_biomass_50km_2019',timeout=40).json()['properties']
    w,_,x0,_,h,y0,*_=meta['proj:transform']
    grid={st:(int((6371007.181*math.radians(r.longitude)*math.cos(math.radians(r.latitude))-x0)//w),
              int((6371007.181*math.radians(r.latitude)-y0)//h)) for st,r in coords.iterrows()}
    env['pixel']=env.station.map(lambda s:str(grid[s]))
    env.to_csv(OUT/'site_year.csv',index=False)
    b=env.dropna(subset=['woody_tonnes'])
    b=b[b.groupby('pixel').land.transform('nunique').eq(1)]
    assert b.groupby(['pixel','year']).woody_tonnes.nunique().max()==1
    b=b.groupby(['pixel','year','land'],as_index=False).agg(woody_tonnes=('woody_tonnes','first'),light_dn=('light_dn','mean'))
    b.to_csv(OUT/'biomass_pixel_year.csv',index=False)
    b['station']=b.pixel
    o=event_panel().merge(env,on=['station','year']).merge(targets[['station','plant']])
    o['unit']=o.station.astype(str)+'|'+o.plant
    o.to_csv(OUT/'flower_onsets.csv',index=False)
    res=fit(env,'summer_ndvi','station')+fit(env,'green_doy','station')+fit(b,'woody_tonnes','pixel')+fit(o,'onset','unit')
    pd.DataFrame(res).to_csv(OUT/'within.csv',index=False)
    eras=[]
    for label,start,end in [('DMSP',1992,2013),('converted_VIIRS',2014,2024)]:
        for frame,metric,group in [(env,'green_doy','station'),(b,'woody_tonnes','pixel')]:
            for row in fit(frame[frame.year.between(start,end)],metric,group):
                eras.append(dict(era=label,**row))
    from scipy.stats import theilslopes
    trends=[]
    for frame,outcome,key in [(env,'green_doy','station'),(b,'woody_tonnes','pixel'),(o,'onset','unit')]:
        for unit,g in frame.dropna(subset=[outcome,'light_dn']).groupby(key):
            if len(g)<3:continue
            trends.append(dict(metric=outcome,unit=unit,land=g.land.iloc[0],first=g.year.min(),last=g.year.max(),
                years=len(g),light_dn_per_year=theilslopes(g.light_dn,g.year).slope,response_per_year=theilslopes(g[outcome],g.year).slope))
    pd.DataFrame(trends).to_csv(OUT/'longitudinal_trends.csv',index=False)
    pd.DataFrame(eras).to_csv(OUT/'sensor_eras.csv',index=False)
    print(pd.DataFrame(res).round(2).to_string(index=False))
if __name__=='__main__':run()
