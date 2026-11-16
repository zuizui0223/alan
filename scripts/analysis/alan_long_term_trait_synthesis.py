"""Finite-cohort precision-weighted synthesis of species slopes with shared-block uncertainty."""
import sys
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests
from alan_global_gbif_screen import model


def traits_for(names, work, results):
    trait = pd.read_csv(results/"regional_light_traits/species_traits.csv",dtype={"sp":str}).set_index("sp").reindex(names.index)
    src = pd.read_csv(work/"direct_trait_replay/direct_species_trait_ledger.csv.gz",low_memory=False)
    src = src[src.quality.isin(["medium","high"]) & src.resolution_status.eq("resolved") & src.multistate.eq(False)]
    axes = src.pivot(index="accepted_species",columns="trait_name",values="normalized_value").reindex(trait.accepted_name).set_axis(names.index)
    def binary(col, yes, no):
        return axes[col].isin(yes).where(axes[col].isin(yes+no)).astype(float)
    out = pd.DataFrame(index=names.index)
    out["annual"] = trait.lifeform.str.contains(r"\bannual\b",case=False).astype(float)
    out["autonomous"] = binary("autonomous_selfing_capacity",["autonomous"],["absent"])
    out["SI"] = binary("self_incompatibility",["SI"],["SC"])
    out["white_tubular"] = (axes.flower_primary_color.eq("white") & axes.floral_form.isin(["tubular","salverform","funnel_trumpet"])).where(axes.flower_primary_color.notna() & axes.floral_form.notna()).astype(float)
    review = trait.transport_review_status.copy()
    review.loc[trait.accepted_name.isin(pd.read_csv(results/"mechanism_exemplars/transport_updates.csv").plant)]="documented"
    out["moth_documented"] = review.eq("documented").where(review.isin(["documented","unconfirmed"])).astype(float)
    out["pollen_limitation"] = trait.pollen_limitation
    return out,trait.wcvp_family


def synthesize(y,scores,covariates,df_min,family,terms,weighted=True):
    x=covariates[terms].dropna().copy()
    x.insert(0,"intercept",1.)
    ids=x.index
    fit={"n":len(ids),"families":family.loc[ids].nunique()}
    if len(terms)==1 and terms[0]!="pollen_limitation":
        fit["yes"]=int(x[terms[0]].sum())
        fit["no"]=len(ids)-fit["yes"]
        if min(fit["yes"],fit["no"])<3:return [dict(fit,term=terms[0],status="sparse")]
    if len(ids)<=x.shape[1]+1 or np.linalg.matrix_rank(x)<x.shape[1]:return []
    z=scores.loc[ids].to_numpy()
    w=1/(z*z).sum(axis=1) if weighted else np.ones(len(ids))
    design=x.to_numpy()
    def projection(use):
        a=design[use]
        return np.linalg.pinv(a.T@(w[use,None]*a))@(a.T*w[use])
    h=projection(np.ones(len(ids),dtype=bool))
    effects=y.loc[ids].to_numpy()
    beta=h@effects
    inf=h@z
    df=df_min.loc[ids].min()
    output=[]
    for j,key in enumerate(x.columns):
        if terms and j==0:continue
        se=np.linalg.norm(inf[j])
        half=t.ppf(.975,df)*se
        omitted=[]
        for fam in family.loc[ids].dropna().unique():
            mask=family.loc[ids].ne(fam).fillna(True).to_numpy(bool)
            if np.linalg.matrix_rank(design[mask])==design.shape[1]:
                omitted.append((projection(mask)@effects[mask])[j])
        output.append(dict(fit,term=key,status="estimated",estimate=beta[j],se=se,low=beta[j]-half,
                           high=beta[j]+half,p=2*t.sf(abs(beta[j]/se),df),df=df,
                           family_omit_min=min(omitted,default=np.nan),family_omit_max=max(omitted,default=np.nan)))
    return output


def main(work,results):
    base=results/"long_term_land_contrasts"
    names=pd.read_csv(base/"all_species.csv.gz",dtype={"sp":str}).dropna(subset=["slope"]).set_index("sp")
    traits,fam=traits_for(names,work,results)
    db=duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    df=db.sql("from trajectories").df()
    effects={key:{} for key in ["pooled","island_minus_mainland","two_land_mean"]}
    for sp,g in df[df.sp.isin(names.index)].groupby("sp"):
        b,se,k,s=model(g)
        np.testing.assert_allclose(b,names.at[sp,"slope"],atol=1e-7)
        effects["pooled"][sp]=(b,k-1,pd.Series(s))
        if pd.isna(names.at[sp,"delta"]):continue
        i=model(g[g.land_class.eq("island")])
        m=model(g[g.land_class.eq("mainland")])
        for key,sign,scale in [("island_minus_mainland",-1,1),("two_land_mean",1,.5)]:
            scores=pd.Series(i[3]).add(sign*pd.Series(m[3]),fill_value=0)*scale
            effects[key][sp]=((i[0]+sign*m[0])*scale,min(i[2],m[2])-1,scores)
        np.testing.assert_allclose(i[0]-m[0],names.at[sp,"delta"],atol=1e-7)
    rows=[]
    for endpoint,values in effects.items():
        y=pd.Series({sp:v[0] for sp,v in values.items()})
        freedom=pd.Series({sp:v[1] for sp,v in values.items()})
        scores=pd.DataFrame({sp:v[2] for sp,v in values.items()}).T.fillna(0)
        covariates=traits.loc[y.index]
        for label,terms,weighted in [("mean",[],True),("equal_species",[],False),
                    *(("trait",[col],True) for col in traits),
                    ("joint",["annual","autonomous","SI"],True)]:
            for row in synthesize(y,scores,covariates,freedom,fam,terms,weighted):
                rows.append(dict(endpoint=endpoint,analysis=label,**row))
    output=pd.DataFrame(rows)
    tests=output.analysis.eq("trait") & output.p.notna()
    output.loc[tests,"q"]=multipletests(output.loc[tests,"p"],method="fdr_bh")[1]
    out=base/"trait_synthesis"
    out.mkdir(exist_ok=True)
    output.to_csv(out/"synthesis.csv",index=False)
    paired=names.dropna(subset=["delta"])
    group=pd.Series(np.select([
        paired.mainland.lt(0)&paired.island.lt(0),
        paired.mainland.lt(0)&paired.island.ge(0),
        paired.mainland.ge(0)&paired.island.lt(0)],
        ["both_negative","mainland_negative","island_negative"],default="both_nonnegative"),
        index=paired.index,name="profile")
    descriptions=[]
    for profile,ids in group.groupby(group).groups.items():
        for trait in traits:
            observed=traits.loc[ids,trait].dropna()
            descriptions.append(dict(profile=profile,species=len(ids),trait=trait,known=len(observed),
                                     positive=observed.sum() if trait!="pollen_limitation" else np.nan,
                                     mean=observed.mean()))
    pd.DataFrame(descriptions).to_csv(out/"profile_composition.csv",index=False)
    complete=traits[["annual","autonomous","SI"]].dropna()
    complete.value_counts().rename("species").to_csv(out/"joint_trait_support.csv")
    print(output[["endpoint","analysis","term","n","yes","estimate","low","high","q"]].round(2).to_string(index=False))


if __name__=="__main__":
    main(*map(Path,sys.argv[1:3]))
