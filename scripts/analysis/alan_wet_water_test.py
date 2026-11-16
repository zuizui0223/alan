"""Held-out prediction for paired wet-plant groups after surface-water adjustment."""
import sys
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import t
from alan_wet_water import stations
from alan_long_term_environment_forest import environmental_months


def observed_brightening(g,out):
    bright=g[g.light>1e-8].copy()
    summary=bright.groupby("state").agg(brightened_coordinates=("unit_id","size"),blocks=("block","nunique"),
        candidate_mean_pp_decade=("candidate_rate","mean"),reference_mean_pp_decade=("reference_rate","mean"),
        paired_mean_pp_decade=("rate_difference","mean"))
    residual=bright.rate_difference-bright.groupby("state").rate_difference.transform("mean")
    scores=residual.groupby([bright.state,bright.block]).sum()
    k=summary.blocks
    se=np.sqrt(scores.pow(2).groupby(level=0).sum()*k/(k-1))/summary.brightened_coordinates
    half=t.ppf(.975,k-1)*se
    summary["paired_low"]=summary.paired_mean_pp_decade-half
    summary["paired_high"]=summary.paired_mean_pp_decade+half
    summary.to_csv(out/"observed_brightening.csv")


def build(work,root,out):
    db=duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    dates=stations(db,root)
    water=pd.read_parquet(work/"gbif_annual/glad_water/survey_water.parquet")
    matched=dates.merge(water,on=["unit_id","y"],validate="1:1")
    assert len(matched)==len(dates)
    matched.groupby("status").size().rename("coordinate_years").to_csv(out/"satellite_coverage.csv")
    water.sort_values(["unit_id","y"]).to_csv(out/"satellite_survey_years.csv.gz",index=False)
    dates=matched[matched.status.eq("ok")].copy()
    db.register("water",dates[["unit_id","y","water_area","water_freq"]])
    db.execute("""create temp table panel as select * from surveys join both_sites using(site) join water using(unit_id,y)
       qualify count(*) over(partition by site)>=5 and max(y) over(partition by site)-min(y) over(partition by site)>=20""")
    db.execute("""create temp table observed as select distinct site,sp,candidate from year_hits join wet using(sp)
                  semi join panel using(site,y)""")
    db.execute("""create temp table sizes as select site,count(*) filter(where candidate) nc,
                     count(*) filter(where not candidate) nr from observed group by site having nc>0 and nr>0""")
    dated=db.sql("""select distinct unit_id,y from panel semi join sizes using(site)""").df()
    d=dates.merge(dated,validate="1:1")
    climate=environmental_months(d.merge(pd.DataFrame({"m":range(1,13)}),how="cross"),work/"gbif_annual/ncep")
    weather=climate.groupby(["unit_id","y"],as_index=False).agg(temp=("air","mean"),rain=("prate","mean"),
                                                              built0=("built0","first"),built=("built_change","first"))
    db.register("weather",weather)
    ss=db.sql("""with hits as (select site,y,sum(hits) filter(where candidate) as hc,
       sum(hits) filter(where not candidate) as hr from year_hits join observed using(site,sp) group by all),
       m as (select p.*,nc,nr,100.*coalesce(hc,0)/(visits*nc) as cr,100.*coalesce(hr,0)/(visits*nr) as rr,
       temp,rain,built0,built,row_number() over(partition by p.site order by p.y) as first_rank,
       row_number() over(partition by p.site order by p.y desc) as last_rank
       from panel p join sizes using(site) join weather using(unit_id,y) left join hits using(site,y))
       select site,unit_id,block,(land_class='island')::int as island,avg(nc) as nc,avg(nr) as nr,min(y) as start,
       max(y)-min(y) as span,count(*) as years,10*regr_slope(cr,y) as candidate_rate,
       10*regr_slope(rr,y) as reference_rate,10*regr_slope(dn,y) as light,arg_min(dn,y) as light0,
       10*regr_slope(water_area,y) as area_trend,10*regr_slope(water_freq,y) as freq_trend,
       avg(water_area) filter(where first_rank<=2) as area0,avg(water_area) filter(where last_rank<=2) as area_end,
       avg(water_freq) filter(where first_rank<=2) as freq0,
       avg(temp) as temp0,avg(rain)*86400 as rain0,10*regr_slope(temp,y) as temp_trend,
       864000*regr_slope(rain,y) as rain_trend,avg(built0) as built0,avg(built) as built,
       10*regr_slope(ln(richness),y) as effort,10*regr_slope(ln(visits),y) as days,
       10*regr_slope(month_fraction_3+month_fraction_4+month_fraction_5,y) as spring,
       10*regr_slope(month_fraction_6+month_fraction_7+month_fraction_8,y) as summer,
       10*regr_slope(month_fraction_9+month_fraction_10+month_fraction_11,y) as autumn
       from m group by all""").df()
    assert ss.site.is_unique
    g=ss.drop(columns="site").groupby(["unit_id","block","island"],as_index=False).mean()
    g["state"]=np.select([g.area0.lt(1),g.area_trend.lt(-2),g.area_trend.gt(2),g.area_end.lt(1)],
                         ["little_initial_open_water","shrinking","expanding","lost_low_cover"],default="maintained")
    g["rate_difference"]=g.candidate_rate-g.reference_rate
    g.to_csv(out/"coordinate_trends.csv.gz",index=False)
    observed_brightening(g,out)
    db.register("states",g[["unit_id","state"]])
    support=db.sql("""select state,candidate,count(distinct sp) as species,count(distinct (sp,unit_id)) as species_coordinates
             from observed join sizes using(site) join geo using(site) join states using(unit_id) group by all""").df()
    support.to_csv(out/"plant_support.csv",index=False)
    db.sql("""select state,sp,candidate,count(distinct unit_id) as coordinates from observed join sizes using(site)
              join geo using(site) join states using(unit_id) group by all""").df().to_csv(out/"species_coverage.csv",index=False)
    species=db.sql("""with stations as (select sp,p.site,unit_id,10*regr_slope(100.*coalesce(hits,0)/visits,y) as rate
        from panel p join sizes using(site) join observed using(site) left join year_hits using(site,y,sp)
        where candidate group by all) select sp,unit_id,avg(rate) as rate from stations group by all""").df()
    species.to_parquet(work/"gbif_annual/glad_water/species_rate.parquet",index=False)
    g.groupby("state").agg(coordinates=("unit_id","size"),blocks=("block","nunique"),island_sites=("island","sum"),
                          area_begin=("area0","mean"),area_end=("area_end","mean"),area_trend=("area_trend","mean")).to_csv(out/"state_counts.csv")
    return g


def predictions(g,scope):
    covs=["temp_trend","rain_trend","built","light0","temp0","rain0","built0","effort","days","spring","summer","autumn",
          "start","span","years","island","nc","nr"]
    water=["area0","freq0","area_trend","freq_trend"]
    g=g.sort_values("unit_id").dropna()
    outputs=[]
    if g.block.nunique()<5:return outputs
    y=g[["candidate_rate","reference_rate","rate_difference"]].to_numpy()
    labels=["candidates","other_wet_species","paired_difference"]
    for repeat in (False,True):
        folds=GroupKFold(5,shuffle=repeat,random_state=391 if repeat else None).split(g,y,g.block)
        errors={m:np.zeros_like(y) for m in ["mean","base","light","water","both"]}
        parts={}
        for i,(tr,te) in enumerate(folds):
            errors["mean"][te]=(y[te]-y[tr].mean(axis=0))**2
            for model,cols in [("base",covs),("light",covs+["light"]),("water",covs+water),("both",covs+water+["light"])]:
                x=g[cols].to_numpy()
                for j in range(3):
                    rf=RandomForestRegressor(n_estimators=64,min_samples_leaf=5,random_state=i+(391 if repeat else 293),n_jobs=1)
                    rf.fit(x[tr],y[tr,j])
                    errors[model][te,j]=(y[te,j]-rf.predict(x[te]))**2
                    if model=="both":
                        light=x[:, -1]
                        edges=np.unique(np.quantile(light[tr],np.linspace(0,1,6)))
                        for a,b in zip(edges[:-1],edges[1:]):
                            ids=te[(light[te]>=a)&((light[te]<b)|((light[te]==b)&(b==edges[-1])))]
                            if len(ids):
                                left=x[ids].copy();right=left.copy()
                                left[:,-1]=a;right[:,-1]=b
                                parts.setdefault(j,[]).extend(zip(g.block.iloc[ids],(rf.predict(right)-rf.predict(left))/(b-a)))
        avg={k:pd.DataFrame(v).groupby(g.block.to_numpy()).mean().mean().to_numpy() for k,v in errors.items()}
        for j,label in enumerate(labels):
            denom=avg["mean"][j]
            outputs.append(dict(scope=scope,outcome=label,split=1+int(repeat),coordinates=len(g),blocks=g.block.nunique(),
                           skill=1-avg["both"][j]/denom,light_gain_before=(avg["base"][j]-avg["light"][j])/denom,
                           light_gain_after=(avg["water"][j]-avg["both"][j])/denom,
                           water_gain_after=(avg["light"][j]-avg["both"][j])/denom,
                           light_direction=pd.DataFrame(parts[j],columns=["block","direction"]).groupby("block").direction.mean().mean()
                           if j in parts else np.nan))
    return outputs


def plot_results(predictions_table,out):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6,3.7))
    colors=["#16786a","#a35232"]
    for i,(label,color) in enumerate(zip(["candidates","other_wet_species"],colors)):
        for j,scope in enumerate(["all_paired","maintained","shrinking"]):
            z=predictions_table.query("outcome==@label and scope==@scope").light_gain_after*100
            if not len(z):continue
            yy=j+(i-.5)*.18
            ax.plot(z,[yy]*len(z),color=color,lw=2)
            ax.scatter(z,[yy]*len(z),color=color,s=30,label=label.replace("_"," ") if j==0 else None)
    ax.set(yticks=[0,1,2],yticklabels=["All paired sites","Maintained open water","Shrinking open water"],
           xlabel="Error reduction from adding light (% of baseline MSE)")
    ax.axvline(0,color=".65",ls="--");ax.legend(frameon=False,fontsize=8)
    fig.tight_layout();fig.savefig(out/"light_after_water.png",dpi=170);plt.close(fig)


def main(work,root):
    out=root/"all_species_forest/wet_water";out.mkdir(exist_ok=True)
    g=build(work,root,out)
    rows=[]
    for scope,df in [("all_paired",g),("maintained",g[g.state.eq("maintained")]),("shrinking",g[g.state.eq("shrinking")])]:
        rows.extend(predictions(df,scope))
    predictions_table=pd.DataFrame(rows)
    predictions_table.to_csv(out/"predictive_tests.csv",index=False)
    plot_results(predictions_table,out)
    spp=pd.read_parquet(work/"gbif_annual/glad_water/species_rate.parquet").merge(g,on="unit_id",validate="m:1")
    roster=pd.read_csv(root/"all_species_forest/decline_traits/moisture_evidence.csv.gz").query("candidate and indicator>7")
    def species_test(sp,group):
        sub=group[group.state.eq("maintained")].copy()
        result=dict(sp=sp,coordinates=len(sub),blocks=sub.block.nunique())
        if len(sub)<10 or sub.block.nunique()<5 or sub.rate.nunique()<2:
            return result
        cols=["temp_trend","rain_trend","built","light0","temp0","rain0","built0","effort","days","spring","summer","autumn",
              "start","span","years","island","area0","freq0","area_trend","freq_trend"]
        for repeat in (False,True):
            testname="repeat" if repeat else "initial"
            kf=GroupKFold(5,shuffle=repeat,random_state=391 if repeat else None)
            errors=np.zeros((len(sub),3))
            signed=[]
            for k,(tr,te) in enumerate(kf.split(sub,groups=sub.block)):
                y=sub.rate.to_numpy()
                errors[te,0]=(y[te]-y[tr].mean())**2
                for j,c in enumerate([cols,cols+["light"]],1):
                    x=sub[c].to_numpy()
                    rf=RandomForestRegressor(n_estimators=32,min_samples_leaf=min(5,max(1,len(tr)//12)),random_state=k+(391 if repeat else 293))
                    rf.fit(x[tr],y[tr]);errors[te,j]=(y[te]-rf.predict(x[te]))**2
                    if j==2:
                        edges=np.unique(np.quantile(x[tr,-1],np.linspace(0,1,6)))
                        for a,b in zip(edges[:-1],edges[1:]):
                            ids=te[(x[te,-1]>=a)&((x[te,-1]<b)|((x[te,-1]==b)&(b==edges[-1])))]
                            if len(ids):
                                l=x[ids].copy();h=l.copy();l[:,-1]=a;h[:,-1]=b
                                signed.extend(zip(sub.block.iloc[ids],(rf.predict(h)-rf.predict(l))/(b-a)))
            m=pd.DataFrame(errors).groupby(sub.block.to_numpy()).mean().mean().to_numpy()
            if m[0]>0:
                result["gain_"+testname]=(m[1]-m[2])/m[0]
                result["skill_"+testname]=1-m[2]/m[0]
                result["slope_"+testname]=pd.DataFrame(signed,columns=["block","effect"]).groupby("block").effect.mean().mean()
        return result
    from joblib import Parallel,delayed
    res=Parallel(n_jobs=5)(delayed(species_test)(sp,group) for sp,group in spp.groupby("sp"))
    table=roster[["sp","species","indicator"]].merge(pd.DataFrame(res),on="sp",how="left",validate="1:1")
    observed=spp[spp.state.eq("maintained") & spp.light.gt(1e-8)].groupby("sp").rate.agg(["size","mean"])
    table.join(observed.rename(columns={"size":"brightening_sites","mean":"brightening_rate"}),on="sp").to_csv(out/"maintained_species.csv",index=False)
    print(g.state.value_counts(),pd.DataFrame(rows).round(3).to_string(index=False))


if __name__=="__main__":main(*map(Path,sys.argv[1:3]))
