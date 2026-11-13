"""Summarize observed animal-plant calendar differences across repeated pairs."""
from alan_decadal_events import *
import json
from scipy.stats import theilslopes

def run():
    d=pd.read_csv(OUT/'species_date_overlap.csv')
    m=pd.read_csv(OUT/'gbif_moth_dates.csv')
    m['year']=pd.to_datetime(m.date).dt.year
    m['day']=pd.to_datetime(m.date).dt.dayofyear
    obs=[]
    for f in C.glob('npn_*_plant.json'):
        obs.extend((r['site_id'],r['species_id'],r['observation_date']) for r in json.loads(f.read_bytes()) if r['phenophase_status']==1)
    o=pd.DataFrame(obs,columns=['station','id','date'])
    o['plant']=o.id.map(names.set_index('species_id').accepted_name)
    o['year']=pd.to_datetime(o.date).dt.year;o['day']=pd.to_datetime(o.date).dt.dayofyear
    flowers=o.groupby(['station','plant','year']).day.agg(set).to_dict()
    rows=[]
    for (st,land,p,s,y),g in m.groupby(['station','land','plant','moth','year']):
        fd=flowers.get((st,p,y));md=set(g.day)
        if fd:rows.append([st,land,p,s,'moth_family_visit',y,len(fd),len(md),len(fd&md),np.median(list(md))-np.median(list(fd))])
    d=pd.concat([d,pd.DataFrame(rows,columns=d.columns)]).drop_duplicates()
    # A bird with both a visit and a pollination record represents one pair.
    d['guild']=d.guild.replace({'pollinates':'bird_flower','visitsFlowersOf':'bird_flower','eats':'bird_food_ripe_fruit_proxy'})
    d=d.drop_duplicates()
    d.to_csv(OUT/'paired_dates.csv',index=False)
    light=pd.read_csv(OUT/'light.csv')
    d=d.merge(light,on=['station','land','year'],validate='m:1')
    fits=[]
    import statsmodels.formula.api as sm
    for guild,g in d.groupby('guild'):
        g=g.copy()
        g['unit']=g.station.astype(str)+'|'+g.plant+'|'+g.animal
        g=g[g.groupby('unit').year.transform('nunique')>=2].copy()
        g['x']=g.light_dn/10
        if g.station.nunique()<3:continue
        mod=sm.ols('animal_minus_plant_median_observed_day ~ x + C(unit) + C(year)',g,eval_env=-1)
        if np.linalg.matrix_rank(mod.exog)<mod.exog.shape[1]:continue
        f=mod.fit(cov_type='cluster',cov_kwds={'groups':g.station},use_t=True)
        fits.append(dict(guild=guild,estimate=f.params.x,lo=f.conf_int().loc['x',0],hi=f.conf_int().loc['x',1],
                         stations=g.station.nunique(),rows=len(g)))
    pd.DataFrame(fits).to_csv(OUT/'paired_date_light.csv',index=False)
    m["doy"]=pd.to_datetime(m.date).dt.dayofyear
    per=m.drop_duplicates(["station","moth","date"]).groupby(["station","land","moth","year"]).doy.agg(["size","min","median","max"]).reset_index()
    per.to_csv(OUT/"moth_calendar_year.csv",index=False)
    ms=[]
    for key,g in per[per["size"]>=3].groupby(["station","land","moth"]):
        if len(g)>=3:ms.append([*key,len(g),g.year.min(),g.year.max(),theilslopes(g["median"],g.year).slope])
    pd.DataFrame(ms,columns=["station","land","moth","years","first","last","observed_median_day_per_year"]).to_csv(OUT/"moth_calendar_trends.csv",index=False)
    slopes=[]
    for keys,g in d.groupby(['station','plant','animal','guild','land']):
        if len(g)>=3:
            slopes.append([*keys,len(g),g.year.min(),g.year.max(),theilslopes(g.animal_minus_plant_median_observed_day,g.year).slope])
    pd.DataFrame(slopes,columns=['station','plant','animal','guild','land','years','first','last','days_per_year']).to_csv(OUT/'paired_date_trends.csv',index=False)
    print(d.groupby('guild').agg(rows=('year','size'),taxa=('animal','nunique'),first=('year','min'),last=('year','max'),same=('same_dates',lambda s:(s>0).sum())))
if __name__=='__main__':run()
