#!/usr/bin/env python3
"""Direct autonomous-selfing portfolio versus prior artificial night light.

This is a global OBSERVATIONAL floral photo-composition analysis, NOT proof
that selfing replaces lost nocturnal pollination on real islands.

True original source: 562 islands, 77,647 0.05-degree photo cell×years,
2,012,925 accepted species×site×years, 1992 to photo-year-minus-one
Li–Cao annual calibrated PCNL; externally published species-direct autonomy
and independent Kendall/Nicholson day/night treatment benefit.

All experimental plant species are REMOVED from the selfing-portfolio response
before computing any night-source co-documentation modifier, to avoid using
one species to define both sides of the interaction.
"""
from __future__ import annotations

import hashlib
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, LIGHT, require, site_outcome_model, fixed_effects_regression, species_model
)
LEDGER_SHA="7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577"
CONTROLS=("loglight","night_presence","light_x_night_presence","log_obs","log_selfing_known")
CLASS_NAMES=("autonomous_selfing_capacity","self_incompatibility")


def direct_plant_traits(ledger,*,high_only=False):
    require(hashlib.sha256(ledger.read_bytes()).hexdigest()==LEDGER_SHA,
            "Independent original source-DIRECT plant trait checksum changed")
    a=pd.read_csv(ledger,dtype=str,keep_default_na=False)
    require((len(a),a.accepted_species.nunique())==(134268,62423),
            "Independent botanical source dimensions changed")
    keep=["high"] if high_only else ["high","medium"]
    z=a[a.trait_name.isin(CLASS_NAMES)
        & a.quality.isin(keep)
        & a.resolution_status.eq("resolved")
        & a.multistate.eq("False")
        & a.normalized_value.ne("")
        & ~a.normalized_value.str.contains("|",regex=False)].copy()
    require(not z.duplicated(["accepted_species","trait_name"]).any(),
            "Conflicted direct trait for accepted botanical species")
    traits=z.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    for k in CLASS_NAMES:
        if k not in traits:traits[k]=pd.NA
    return traits[list(CLASS_NAMES)]


def portfolio_frame(d,sites,experiment,trait_table,trait):
    """Retain measured yes/no, not arbitrary 'selfing unknown' as no."""
    require(trait in CLASS_NAMES,"Unknown trait")
    studied=set(experiment.species)
    experimental_night=set(experiment.loc[experiment.night_stronger,"species"])
    experimental_day=set(experiment.loc[experiment.day_stronger,"species"])

    p=d[KEY+["species","genus"]].copy()
    # Separate independent source labels for night-better and daytime-better.
    p["night"] = p.species.isin(experimental_night)
    p["day"] = p.species.isin(experimental_day)
    x=p[KEY+["night","day"]].groupby(KEY,observed=True).sum().reset_index()
    q=sites.merge(x,on=KEY,validate="one_to_one")
    q["night_presence"]=q.night.gt(0).astype(int)
    q["day_presence"]=q.day.gt(0).astype(int)

    p=p.loc[~p.species.isin(studied)].copy()  # independent botanicals
    p["source_direct_trait"]=p.species.map(trait_table[trait])
    yes,no=(("autonomous","absent") if trait=="autonomous_selfing_capacity"
            else ("SC","SI"))
    p["trait_yes"]=p.source_direct_trait.eq(yes).fillna(False)
    p["trait_no"]=p.source_direct_trait.eq(no).fillna(False)
    p["trait_known"]=p.trait_yes|p.trait_no
    g=p.groupby(KEY,observed=True)[["trait_yes","trait_no","trait_known"]].sum().reset_index()
    q=q.merge(g,on=KEY,how="left",validate="one_to_one")
    for k in ("trait_yes","trait_no","trait_known"):
        q[k]=q[k].fillna(0).astype(int)
    require((q.trait_yes+q.trait_no==q.trait_known).all(),
            "Direct trait denominator includes unclassified species")
    require((q.trait_known<=q.all_accepted_recorded_species).all(),
            "Direct-trait coverage exceeds true source flora")
    q["study_experimental_source_outcome_excluded"]=True
    return q,p


def compatible(s,*,need_night_variation=False):
    """Match same original island x YEAR and genuinely distinct light sites."""
    q=s.loc[s.trait_known.ge(2)].copy()
    original_n=int(len(q));islands_orig=int(q.island.nunique())
    meta=dict(
        known_eligible_source_sites_before_light=original_n,
        islands_with_known_eligibility=islands_orig,
        source_sites_night_study_positive_before_light=int(q.night_presence.sum()),
        source_sites_day_study_positive_before_light=int(q.day_presence.sum()))
    if q.empty:
        return q,meta
    g=q.groupby(["island","year"],observed=True).agg(
        n_sites=(LIGHT,"size"),light_min=(LIGHT,"min"),light_max=(LIGHT,"max"),
        night_classes=("night_presence","nunique"))
    good=(g.n_sites.ge(2) & (g.light_max-g.light_min).ge(.02))
    if need_night_variation:
        good=good & g.night_classes.eq(2)
    valid=g.loc[good]
    q=q.merge(valid[[]].reset_index(),on=["island","year"],validate="many_to_one")
    meta.update({
        "eligible_island_years":int(len(valid)),
        "eligible_islands":int(q.island.nunique()),
        "eligible_photo_site_years":int(len(q)),
        "positive_night_study_source_photo_sites":int(q.night_presence.sum()),
        "negative_night_study_documentation_sites":int((1-q.night_presence).sum()),
        "signed_light_span_min_per_island_year":.02,
    })
    return q,meta


def weighted_multivariate_year_island_FE(s,response,*,weight="island_equal"):
    """Weighted island-year FE, island-cluster CR1, t with G−1 df.

    Returns all controls, including night*light. No GloBI absence treated
    as absence of a true nighttime pollinator.
    """
    if len(s)<500 or s.island.nunique()<30 or s.night_presence.nunique()!=2:
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_ISLAND_SOURCE_SUPPORT"}
    g,group=pd.factorize(s.island.astype(str).str.cat(
        s.year.astype(str),sep="|"),sort=True)
    islands,levels=pd.factorize(s.island,sort=True)
    w=(1/s.groupby("island",observed=True).island.transform("size").to_numpy(float)
       if weight=="island_equal" else np.ones(len(s)))
    y=np.asarray(response,float)
    L=s[LIGHT].to_numpy(float)
    N=s.night_presence.to_numpy(float)
    L=L-np.median(L)  # centering only; contrasts remain original L unit
    X=np.column_stack([L,N,L*N,s.log_obs.to_numpy(float),
                      np.log1p(s.trait_known.to_numpy(float))])
    names=list(CONTROLS)
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("Non-finite original source covariates")
    wsum=np.bincount(g,weights=w,minlength=len(group))
    sw=np.sqrt(w)
    for j in range(len(names)):
        mu=np.bincount(g,weights=w*X[:,j],minlength=len(group))/wsum
        X[:,j]=(X[:,j]-mu[g])*sw
    ymu=np.bincount(g,weights=w*y,minlength=len(group))/wsum
    yy=(y-ymu[g])*sw
    rank=int(np.linalg.matrix_rank(X))
    if rank!=X.shape[1]:
        return {"status":"NOT_ESTIMABLE_SOURCE_FIXED_EFFECT_RANK",
                "regressor_rank":rank,"expected_rank":int(X.shape[1])}
    bread=np.linalg.inv(X.T@X)
    coef=bread@X.T@yy
    resid=yy-X@coef
    G=len(levels)
    scores=np.zeros((G,len(names)))
    np.add.at(scores,islands,X*resid[:,None])
    degfree=len(s)-len(names)-len(group)
    if degfree<=0:return {"status":"NOT_ESTIMABLE_FIXED_EFFECT_DEGREES"}
    cov=(G/(G-1)) * ((len(s)-1)/degfree) * bread@scores.T@scores@bread
    se=np.sqrt(np.maximum(0.,np.diag(cov)))
    crit=student_t.ppf(.975,G-1)
    res={k:{"beta_pp":float(coef[i]*100),
            "95CI_pp":[float(100*(coef[i]-crit*se[i])),
                       float(100*(coef[i]+crit*se[i]))],
            "p_nominal":float(2*student_t.sf(abs(coef[i]/se[i]),G-1))
                       if se[i]>0 else None}
         for i,k in enumerate(names)}
    return {"status":"EXPLORATORY_PREVIOUS_LIGHT_X_NIGHT_EXPERIMENT_STUDY_DOCUMENTATION",
            "n_islands":G,"n_source_photo_site_years":len(s),
            "n_island_x_year_fixed_effects":len(group),
            "weighting":weight,"coefficients":res,
            "original_experiment_species_excluded_from_direct_selfing_outcome":True}


def run_selfing_portfolio(d,sites,experiment,ledger):
    res={}
    for quality in ("direct_high_medium","direct_high_only"):
        table=direct_plant_traits(ledger,high_only=(quality=="direct_high_only"))
        q,p=portfolio_frame(d,sites,experiment,table,"autonomous_selfing_capacity")
        meta={
            "external_direct_trait_species":int(table.autonomous_selfing_capacity.isin(
                ["autonomous","absent"]).sum()),
            "photo_source_flora_species_direct_known_after_study_taxon_exclusion":
                int(p.loc[p.trait_known,"species"].nunique()),
            "photo_source_flora_species_direct_yes_after_study_taxon_exclusion":
                int(p.loc[p.trait_yes,"species"].nunique()),
            "photo_source_flora_species_direct_no_after_study_taxon_exclusion":
                int(p.loc[p.trait_no,"species"].nunique()),
            "source_selfing_known_species_site_years":int(p.trait_known.sum()),
            "source_selfing_yes_species_site_years":int(p.trait_yes.sum()),
        }
        # Include documentation-selection sources; these are NOT biological
        # not-selfing control species.
        diag={}
        for name in ("trait_known","trait_yes","trait_no"):
            diag[name]=site_outcome_model(q,q[name].to_numpy(float))
        study_all,all_support=compatible(q,need_night_variation=False)
        if (len(study_all)<500 or study_all.island.nunique()<30):
            proportion={"status":"NOT_ESTIMABLE_KNOWN_SELFING_GEOGRAPHIC_SUPPORT",
                        "source_support":all_support}
        else:
            y=study_all.trait_yes.to_numpy(float)/study_all.trait_known.to_numpy(float)
            w=1/study_all.groupby("island").island.transform("size").to_numpy(float)
            iy=study_all.island.astype(str).str.cat(study_all.year.astype(str),sep="|")
            fit=fixed_effects_regression(y,study_all[LIGHT].to_numpy(float),
                     study_all.log_obs.to_numpy(float),w,study_all.island.to_numpy(),
                     [iy.to_numpy()])
            proportion={"status":"EXPLORATORY_DIRECT_SELFING_KNOWN_TAXA_ONLY",
                        "fit":fit,"source_support":all_support,
                        "do_not_infer_selfing_of_species_without_DIRECT_TRAIT":True}
        cand,cond_support=compatible(q,need_night_variation=True)
        if (len(cand)<500 or cand.island.nunique()<30 or
            min(int(cand.night_presence.sum()),int((1-cand.night_presence).sum()))<100):
            interaction={"status":"NOT_ESTIMABLE_NIGHT_STUDY_SOURCE_STRATA_SUPPORT",
                         "source_support":cond_support}
        else:
            yy=cand.trait_yes.to_numpy(float)/cand.trait_known.to_numpy(float)
            interaction={
                "island_equal":weighted_multivariate_year_island_FE(cand,yy),
                "site_year_equal":weighted_multivariate_year_island_FE(
                    cand,yy,weight="site_equal"),
                "source_support":cond_support,
                "modifier_is_OBSERVED_EXTERNAL_NIGHT_STUDY_TAXON_NOT_NOCTURNAL_POLLINATOR_PRESENCE":True
            }
        # Secondary conceptual negative control with confirmed self-compatibility.
        qc,pc=portfolio_frame(d,sites,experiment,table,"self_incompatibility")
        sc,scmeta=compatible(qc,need_night_variation=True)
        if (len(sc)<500 or sc.island.nunique()<30 or
            min(sc.night_presence.sum(),len(sc)-sc.night_presence.sum())<100):
            scfit={"status":"NOT_ESTIMABLE_SC_SI_STRATA_SUPPORT",
                   "source_support":scmeta}
        else:
            val=sc.trait_yes.to_numpy(float)/sc.trait_known.to_numpy(float)
            scfit={"result":weighted_multivariate_year_island_FE(sc,val),
                   "source_support":scmeta}
        res[quality]={
           "source_species_direct_selfing_coverage":meta,
           "source_ascertainment_controls":diag,
           "observed_autonomous_selfing_ratio_among_directly_known":proportion,
           "external_night_reproduction_species_presence_moderation":interaction,
           "independent_sc_si_documentation_control":scfit
        }
    return res
