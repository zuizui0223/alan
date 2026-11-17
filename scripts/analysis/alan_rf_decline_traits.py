"""Trait contrasts for the locked 283 light-associated declines."""
import sys
from pathlib import Path
from zipfile import ZipFile
import hashlib
import duckdb
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def annotate(work,results,frame):
    eiv=pd.read_csv(work/"gbif_annual/EIV_Europe.csv",index_col=0,encoding="latin1").apply(pd.to_numeric,errors="coerce")
    with ZipFile(work/"v011/wcvp.zip") as z:
        w=pd.read_csv(z.open("wcvp_names.csv"),sep="|",dtype=str,usecols=["taxon_name","plant_name_id","accepted_plant_name_id","taxon_status","taxon_rank","family","lifeform_description"])
    w=w[w.taxon_rank.eq("Species")]
    accepted=w[w.taxon_status.eq("Accepted")].set_index("plant_name_id")
    gl=pd.read_csv(work/"pollen/gl.csv",encoding="latin1")
    gl["name"]=gl.Species_accepted_names.str.replace("_"," ",regex=False).str.strip()
    wanted=set(frame.species)|set(eiv.index)|set(gl.name.dropna())
    aliases=w[w.taxon_name.isin(wanted)].copy()
    aliases["id"]=aliases.accepted_plant_name_id.where(~aliases.taxon_status.eq("Accepted"),aliases.plant_name_id)
    aliases["valid"]=aliases.taxon_status.isin(["Accepted","Synonym","Heterotypic Synonym","Homotypic Synonym"]) & aliases.id.isin(accepted.index)
    groups=aliases.groupby("taxon_name")
    name_id=groups.id.first().where(groups.valid.all() & groups.id.nunique().eq(1)).dropna()
    ids=frame.species.map(name_id)
    name=ids.map(accepted.taxon_name)
    family=ids.map(accepted.family)
    life=ids.map(accepted.lifeform_description)
    t=pd.DataFrame(index=frame.index)
    for col,pattern in {"annual":r"\bannual\b","woody":"tree|shrub|liana","geophyte":"geophyte","aquatic_lifeform":"hydrophyte|helophyte"}.items():
        t[col]=life.str.contains(pattern,case=False).astype(float)
    eiv.index=eiv.index.map(name_id)
    E=eiv.groupby(level=0)
    clean=E.first().where(E.nunique()<=1)
    for col in clean:t[col]=ids.map(clean[col])
    axes=pd.read_csv(work/"direct_trait_replay/direct_species_trait_ledger.csv.gz",low_memory=False)
    axes=axes[axes.quality.isin(["medium","high"])&axes.resolution_status.eq("resolved")&axes.multistate.eq(False)]
    a=axes.pivot(index="accepted_species",columns="trait_name",values="normalized_value").reindex(name).set_axis(frame.index)
    for col,yes,no in [("autonomous_selfing_capacity","autonomous","absent"),("self_incompatibility","SI","SC")]:
        t[col]=a[col].eq(yes).where(a[col].isin([yes,no])).astype(float)
    t["predominantly_selfing"]=a.mating_system.isin(["predominantly_selfing","obligate_selfing"]).where(a.mating_system.notna()).astype(float)
    t["white_tubular"]=(a.flower_primary_color.eq("white") & a.floral_form.isin(["tubular","salverform","funnel_trumpet"])).where(a.flower_primary_color.notna() & a.floral_form.notna()).astype(float)
    ledger=pd.read_csv(results/"species_interpretation/shared_544.csv").set_index("plant")
    review=name.map(ledger.transport_review_status)
    review.loc[name.isin(pd.read_csv(results/"mechanism_exemplars/transport_updates.csv").plant)]="documented"
    t["moth_transport_reviewed"]=review.eq("documented").where(review.isin(["documented","unconfirmed"])).astype(float)
    t["night_exclusion_loss"]=name.map(ledger.day_only_relative_loss)
    t["pollen_limitation"]=ids.map(gl.groupby(gl.name.map(name_id)).PL_Effect_Size.median())
    flags=pd.read_csv(next((work/"v011").glob("*flags/*csv.gz")))
    cross=pd.read_csv(work/"v011/moth_name_crosswalk.csv")
    positive=flags.loc[flags.non_inat_visits_moth_species_insect.eq(1),"plant_name"]
    moth=cross.loc[cross.resolution.eq("resolved")&cross.name.isin(positive),"accepted_name"]
    t["moth_visit_documented"]=name.isin(moth).where(name.notna()).astype(float)
    return t,family,name


def compare(frame,values,keys,rng):
    n=0; matched_mean=0.; ref_mean=0.; sims=np.zeros(3999); strata=0
    for _,g in frame.assign(v=values).dropna(subset=["v",*keys]).groupby(keys,observed=True):
        a=g.loc[g.candidate,"v"].to_numpy(); b=g.loc[~g.candidate,"v"].to_numpy()
        if not len(a) or not len(b):continue
        k=len(a); n+=k; strata+=1
        matched_mean+=a.sum(); ref_mean+=k*b.mean()
        v=np.r_[a,b]; size=len(v)
        # Preserve the case count in each geographic/support stratum.
        for begin in range(0,len(sims),200):
            chunk=sims[begin:begin+200]
            idx=np.argpartition(rng.random((len(chunk),size)),k-1,axis=1)[:,:k]
            s=v[idx].sum(axis=1)
            chunk+=s-k*(v.sum()-s)/(size-k)
    row=dict(known=int(values[frame.candidate].notna().sum()),matched=n,strata=strata)
    if n:
        diff=matched_mean-ref_mean
        row.update(candidate_mean=matched_mean/n,reference_mean=ref_mean/n,difference=diff/n,
                   p=(1+np.count_nonzero(abs(sims)>=abs(diff)-1e-9))/(len(sims)+1))
    return row


def run(work,results):
    root=results/"all_species_forest";out=root/"decline_traits";out.mkdir(exist_ok=True)
    frame=pd.read_csv(root/"all_species.csv.gz").query("status=='evaluated'").set_index("sp")
    cases=pd.read_csv(root/"decline_candidates_broad.csv").sp
    frame["candidate"]=frame.index.isin(cases)
    assert frame.candidate.sum()==283 and len(frame)==3901
    db=duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    xy=db.sql("select distinct unit_id,longitude,latitude from geo").df()
    raw=pd.read_parquet(work/"gbif_annual/all_species_rf_sites.parquet",columns=["sp","unit_id","island"]).merge(xy,on="unit_id",validate="m:1")
    raw["region"]=(np.floor((raw.longitude+180)/5)*100+np.floor((raw.latitude+90)/5)).astype(int)
    counts=raw.groupby(["sp","region"]).size().reset_index(name="n").sort_values(["sp","n","region"],ascending=[True,False,True])
    frame["region"]=counts.drop_duplicates("sp").set_index("sp").region
    frame["island_band"]=pd.cut(raw.groupby("sp").island.mean(),[-.01,.25,.75,1],labels=False)
    frame["site_band"]=pd.cut(frame.sites,[0,10,30,100,300,1000,np.inf],labels=False)
    t,frame["family"],name=annotate(work,results,frame)
    rows=[];rng=np.random.default_rng(27183)
    for group,subset,keys in [("all",frame,["region","island_band","site_band"]),
          ("other_declines",frame[frame.observed_at_brightening.lt(0)],["region","island_band","site_band"]),
          ("within_family",frame,["region","island_band","site_band","family"])]:
        for col in t:
            row=compare(subset,t.loc[subset.index,col],keys,rng)
            rows.append(dict(comparison=group,trait=col,**row))
    d=pd.DataFrame(rows)
    for _,grp in d.dropna(subset=["p"]).groupby("comparison"):
        d.loc[grp.index,"q"]=multipletests(grp.p,method="fdr_bh")[1]
    d.to_csv(out/"traits.csv",index=False)
    # Diagnose whether the strongest habitat contrast depends on taxonomy or one region.
    extra=[]
    exact=pd.read_csv(work/"gbif_annual/EIV_Europe.csv",index_col=0,encoding="latin1").MOIST
    exact=pd.to_numeric(exact,errors="coerce").groupby(level=0).median()
    direct=frame.species.map(exact)
    frame["light_band"]=pd.qcut(pd.read_parquet(work/"gbif_annual/all_species_rf_sites.parquet",
                        columns=["sp","light0"]).groupby("sp").light0.mean().reindex(frame.index),
                        3,labels=False,duplicates="drop")
    for label,sub,moist,keys in [
        ("exact_name",frame,direct,["region","island_band","site_band"]),
        ("baseline_light",frame,t.MOIST,["region","island_band","site_band","light_band"]),
        ("island_heavy",frame[frame.island_band.eq(2)],t.MOIST,["region","site_band"]),
        ("mainland_heavy",frame[frame.island_band.eq(0)],t.MOIST,["region","site_band"]),
        ("mixed",frame[frame.island_band.eq(1)],t.MOIST,["region","site_band"])]:
        extra.append(dict(check=label,**compare(sub,moist.loc[sub.index],keys,rng)))
    pd.DataFrame(extra).to_csv(out/"moisture_checks.csv",index=False)
    pd.DataFrame(dict(species=frame.species,candidate=frame.candidate,indicator=t.MOIST,
                      site_region=frame.region,island_band=frame.island_band,
                      site_band=frame.site_band)).to_csv(out/"moisture_evidence.csv.gz")
    strata=[]
    bins=[-np.inf,3,5,7,9,np.inf]
    observed=np.zeros(5); expected=np.zeros(5)
    for (region,_,_),g in frame.assign(v=t.MOIST).dropna(subset=["v"]).groupby(["region","island_band","site_band"]):
        a=g.loc[g.candidate,"v"];b=g.loc[~g.candidate,"v"]
        if a.empty or b.empty:continue
        strata.append((region,len(a),a.sum()-len(a)*b.mean()))
        observed+=pd.cut(a,bins).value_counts(sort=False).to_numpy()
        expected+=len(a)*pd.cut(b,bins).value_counts(normalize=True,sort=False).to_numpy()
    block=pd.DataFrame(strata,columns=["region","n","contrast_sum"]).groupby("region").sum()
    block["difference_without_region"]=(block.contrast_sum.sum()-block.contrast_sum)/(block.n.sum()-block.n)
    block.to_csv(out/"moisture_region_influence.csv")
    np.testing.assert_allclose(observed.sum(),expected.sum())
    mat=pd.DataFrame({"score":["<=3","(3,5]","(5,7]","(7,9]",">9"],"candidate_percent":100*observed/observed.sum(),
                      "reference_percent":100*expected/expected.sum()})
    mat.to_csv(out/"moisture_distribution.csv",index=False)
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6.5,3.8))
    p=np.arange(5)
    ax.bar(p-.2,mat.candidate_percent,.38,color="#19655c",label="Decline candidates")
    ax.bar(p+.2,mat.reference_percent,.38,color="#acb5ba",label="Matched reference")
    ax.set(xticks=p,xticklabels=mat.score,xlabel="Habitat moisture indicator (higher = wetter)",ylabel="Species (%)")
    ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(out/"moisture.png",dpi=180);plt.close(fig)
    # Direct source trait classifications remain access-controlled; publish only aggregates.
    frame.loc[frame.candidate,["species","family","sites","blocks","region","island_band","site_band"]].assign(accepted_name=name).to_csv(out/"cohort.csv")
    print(d.query("comparison=='all'").round(3).to_string(index=False))
    print("eiv sha256",hashlib.sha256((work/"gbif_annual/EIV_Europe.csv").read_bytes()).hexdigest())


if __name__=="__main__":run(*map(Path,sys.argv[1:3]))
