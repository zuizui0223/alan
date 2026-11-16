"""Compare annual plant recording rates on identical surveys with both reproductive groups."""
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import t
from scipy.optimize import minimize_scalar
from scipy import sparse


def species_model(db, df):
    """Random species slopes for the annual difference of two roster means."""
    memberships = db.sql("select * from typed").df()
    rows = []
    for land, g in df.groupby("land_class"):
        g = g.reset_index(drop=True)
        strata, keys = pd.factorize(pd.MultiIndex.from_frame(g[["site","m"]]))
        members = memberships.merge(g[["site","m"]].assign(k=strata).drop_duplicates(),
                                    on=["site","m"],validate="m:1")
        sp, taxa = pd.factorize(members.sp)
        n = members.groupby(["k","normalized_value"]).sp.transform("size")
        weights = np.where(members.normalized_value.eq("autonomous"),1,-1)/n
        roster = sparse.csr_matrix((weights,(members.k,sp)),shape=(len(keys),len(taxa)))
        dn = (g.dn-g.groupby(["site","m"]).dn.transform("mean")).to_numpy()/10
        z = roster[strata].multiply(dn[:,None]).tocsr()
        # Project out station-month and year effects using sparse normal equations.
        fy = pd.get_dummies(g.y,dtype=float).to_numpy()
        year = fy-pd.DataFrame(fy).groupby(strata).transform("mean").to_numpy()
        effort = (g.effort-g.groupby(["site","m"]).effort.transform("mean")).to_numpy()
        x = np.column_stack([year,effort,dn])
        y = g.aut-g.absent
        y = (y-y.groupby(strata).transform("mean")).to_numpy()
        xx = x.T@x
        inv = np.linalg.pinv(xx,hermitian=True)
        zx = z.T@x
        zy = z.T@y
        zz = (z.T@z).toarray()
        eig, vec = np.linalg.eigh(zz-zx@inv@zx.T)
        keep = eig > max(1e-9,eig.max()*1e-10)
        eig = eig[keep]
        proj = vec[:,keep].T@(zy-zx@inv@(x.T@y))
        rss = np.sum((y-x@(inv@(x.T@y)))**2)
        dof = len(y)-len(keys)-np.linalg.matrix_rank(x)
        def profile(ratio):
            adj = rss-np.sum(proj**2*ratio/(1+ratio*eig))
            return dof*np.log(adj/dof)+np.log1p(ratio*eig).sum()
        opt = minimize_scalar(lambda log_r:profile(np.exp(log_r)),bounds=(-18,18),method="bounded")
        ratio = np.exp(opt.x) if opt.fun < profile(0) else 0
        vals, vect = np.linalg.eigh(zz)
        def vinv(v):
            return v-z@(vect@((ratio/(1+ratio*np.maximum(vals,0)))[:,None]*(vect.T@(z.T@v))))
        vi_x = vinv(x)
        h = np.linalg.pinv(x.T@vi_x,hermitian=True)
        beta = h@(vi_x.T@y)
        # Independent fixed-effect solution of the penalized mixed-model equations.
        if ratio > 0:
            from scipy.linalg import solve
            shrink = solve(zz + np.eye(len(taxa))/ratio, np.column_stack([zx,zy]),assume_a="pos")
            check = np.linalg.pinv(xx-zx.T@shrink[:,:-1],hermitian=True)@(x.T@y-zx.T@shrink[:,-1])
            np.testing.assert_allclose(check[-1],beta[-1],rtol=1e-6,atol=1e-8)
        block = pd.factorize(g.block)[0]
        influence = (vi_x@h[:,-1])*(y-x@beta)
        score = np.bincount(block,weights=influence)
        sigma2 = (rss-np.sum(proj**2*ratio/(1+ratio*eig)))/dof
        rows.append(dict(land=land,species=len(taxa),random_slope_sd=np.sqrt(ratio*sigma2),
                         **interval(beta[-1],score)))
    return pd.DataFrame(rows)


def residualize(matrix, codes):
    a = matrix.copy()
    for _ in range(3000):
        change = 0
        for code in codes:
            count = np.bincount(code)
            for col in a.T:
                mean = np.bincount(code, weights=col) / count
                col -= mean[code]
                change = max(change, abs(mean).max())
        if change < 1e-9:
            return a
    raise RuntimeError("Fixed effects failed to converge")


def fit(g):
    strata = pd.factorize(pd.MultiIndex.from_frame(g[["site","m"]]))[0]
    years = pd.factorize(g.y)[0]
    block, ids = pd.factorize(g.block)
    a = residualize(g[["dn","effort","aut","absent"]].to_numpy(float), [strata,years])
    a[:,0] /= 10
    light = a[:,0] - a[:,1] * np.linalg.lstsq(a[:,1:2],a[:,0],rcond=None)[0][0]
    info = light @ light
    if info < 1e-8:
        return None
    b = light @ a[:,2:] / info
    error = a[:,2:] - a[:,:2] @ np.linalg.lstsq(a[:,:2],a[:,2:],rcond=None)[0]
    scores = pd.DataFrame(error*light[:,None]/info).groupby(block).sum().set_axis(ids)
    return b, scores


def interval(b,scores):
    k=len(scores)
    se=np.sqrt((scores**2).sum()*k/(k-1))
    half=t.ppf(.975,k-1)*se
    return dict(slope=b,se=se,low=b-half,high=b+half,blocks=k)


def main(work):
    root=Path(__file__).resolve().parents[2]/"results/regional_light_traits"
    out=root/"annual_autonomy"
    out.mkdir(exist_ok=True)
    names=pd.read_csv(root/"species_traits.csv")
    a=pd.read_csv(work/"direct_trait_replay/direct_species_trait_ledger.csv.gz",low_memory=False)
    a=a[a.trait_name.eq("autonomous_selfing_capacity") & a.quality.isin(["medium","high"])
        & a.resolution_status.eq("resolved") & a.multistate.eq(False)
        & a.normalized_value.isin(["autonomous","absent"])]
    roster=names.merge(a[["accepted_species","normalized_value"]],left_on="accepted_name",right_on="accepted_species",validate="m:1")
    db=duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    screened = pd.read_csv(root/"all_regions.csv.gz").query("blocks >= 3")
    traits = screened[["sp","region"]].merge(roster[["sp","normalized_value"]],validate="m:1")
    db.register("regional_traits",traits)
    db.execute("""create temp table locations as select site,case when land_class='island' then 'I_'||land_id else
        'M_'||cast(floor((longitude+180)/5) as int)||'_'||cast(floor((latitude+90)/5) as int) end as region from geo""")
    db.execute("""create temp table typed as select site,m,sp,normalized_value
        from risk join locations using(site) join regional_traits using(sp,region)""")
    db.execute("""create temp table rosters as select site,m,
        count(*) filter(where normalized_value='autonomous') as na,
        count(*) filter(where normalized_value='absent') as nb
        from typed group by all having na>0 and nb>0""")
    df=db.sql("""with tally as (select site,m,y,
        sum(hits) filter(where normalized_value='autonomous') as ha,
        sum(hits) filter(where normalized_value='absent') as hb
        from hits join typed using(site,m,sp) group by all),
        varying as (select site,m from supported group by all having max(dn)>min(dn))
        select s.site,s.m,s.y,s.dn,s.block,s.land_class,s.unit_id,s.visits,na,nb,
        ln(s.richness) as effort,100.*coalesce(ha,0)/(na*s.visits) as aut,
        100.*coalesce(hb,0)/(nb*s.visits) as absent,
        case when s.land_class='island' then 'I_'||s.land_id else
        'M_'||cast(floor((longitude+180)/5) as int)||'_'||cast(floor((latitude+90)/5) as int) end as region
        from supported s join rosters using(site,m) join varying using(site,m)
        join geo g on g.site=s.site left join tally on tally.site=s.site and tally.m=s.m and tally.y=s.y""").df()
    assert not df.duplicated(["site","m","y"]).any() and df[["aut","absent"]].stack().between(0,100).all()
    species_model(db,df).to_csv(out/"hierarchical_species.csv",index=False)
    db.register("cells",df[["site","m","land_class"]].drop_duplicates())
    counts=db.sql("""select land_class,normalized_value,count(distinct sp) as species from cells
        join typed using(site,m) group by all""").df()
    cov=df.groupby("land_class").agg(month_years=("y","size"),survey_days=("visits","sum"),coordinates=("unit_id","nunique"),
        years=("y","nunique"),first_year=("y","min"),last_year=("y","max"),regions=("region","nunique"),blocks=("block","nunique"))
    cov.join(counts.pivot(index="land_class",columns="normalized_value",values="species")).to_csv(out/"coverage.csv")
    rows=[]
    scores={}
    for land,g in df.groupby("land_class"):
        b,s=fit(g)
        for i,label in enumerate(["autonomous","absent"]):
            rows.append(dict(land=land,contrast=label,**interval(b[i],s[i])))
        delta=b[0]-b[1]
        sc=s[0]-s[1]
        rows.append(dict(land=land,contrast="autonomous_minus_absent",**interval(delta,sc)))
        scores[land]=(delta,sc)
    bi,si=scores["island"]
    bm,sm=scores["mainland"]
    rows.append(dict(land="island_minus_mainland",contrast="autonomous_minus_absent",
                     **interval(bi-bm,si.subtract(sm,fill_value=0))))
    pd.DataFrame(rows).to_csv(out/"annual_contrasts.csv",index=False)
    regional=[]
    for (land,region),g in df.groupby(["land_class","region"]):
        row=dict(land=land,region=region,coordinates=g.unit_id.nunique(),blocks=g.block.nunique(),cells=len(g))
        if row["blocks"]>=3:
            f=fit(g)
            if f:
                b,s=f
                row.update(interval(b[0]-b[1],s[0]-s[1]))
        regional.append(row)
    regional=pd.DataFrame(regional)
    regional.to_csv(out/"regional.csv",index=False)
    # Random regional slopes, with uncertainty from the first-stage annual models.
    pooled=[]
    for land,g in regional.dropna(subset=["se"]).groupby("land"):
        g=g[g.se>1e-7]
        if len(g)<3: continue
        b=g.slope.to_numpy()
        v=g.se.to_numpy()**2
        def objective(tau2):
            w=1/(v+tau2)
            return (np.log(v+tau2).sum()+np.log(w.sum())+np.sum(w*(b-np.average(b,weights=w))**2))/2
        opt=minimize_scalar(objective,bounds=(0,max(1,np.ptp(b)**2)),method="bounded")
        tau2=opt.x if opt.fun<objective(0) else 0
        w=1/(v+tau2)
        mean=np.average(b,weights=w)
        se=np.sqrt(max(1,np.sum(w*(b-mean)**2)/(len(g)-1))/w.sum())
        half=t.ppf(.975,len(g)-1)*se
        pooled.append(dict(land=land,regions=len(g),slope=mean,low=mean-half,high=mean+half,tau=np.sqrt(tau2)))
    pd.DataFrame(pooled).to_csv(out/"hierarchical_regions.csv",index=False)


if __name__=="__main__":
    main(Path(sys.argv[1]))
