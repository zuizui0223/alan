#!/usr/bin/env python3
"""ALAN: decompose observed island-size × nightlight into photographed taxon mix
and taxon-relative ORIGINAL photo-ROI white components.

Same 2019/2021/2023/2024/2025 ROI9-source/GSHHG/WCVP/VIIRS/WorldPop inputs,
same 44 strict or 95 broad actually photographed multi-site islands and same
island-equal/cell-equal regression weights as the previously signed area and
mainland-isolation study. No additional images or synthetic ecological labels.

The full-source taxon reference is an *internal photo archive mean*, not an
independent botanically verified species trait: full-sample decomposition is
an algebraic partition, NOT causal species turnover. A separate holdout
learns source taxon photo colour in 2019/21/23 then tests 2024/25 without
using those target photos in the reference. The holdout has a DIFFERENT
observed-species reference-coverage subset, explicitly documented.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from scripts.analysis import alan_island_area_light_20261108 as area
from scripts.analysis import alan_isolation_area_light_20261108 as isolation
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis import mois_island_ALAN_taxon_turnover_measured_photo_20261107 as taxa

YEARS=(2019,2021,2023,2024,2025)
PRIOR=(2019,2021,2023)
HOLDOUT=(2024,2025)
COHORT_N={"strict_5p3r3s":7982,"measured_3p2r2s":20543}
SOURCE_ISLAND_N={"strict_5p3r3s":44,"measured_3p2r2s":95}
SOURCE_CELL_N={"strict_5p3r3s":650,"measured_3p2r2s":1884}
BASE_COVARS="log_pop + bio5 + log_rain + elevation + log_attempts"
TERMS=("log_light:log10_area_c","log_light:log10_distance_c")
ORIGINAL_BETAS={
  "strict_5p3r3s":{
     "cell_equal":(2.216703936269394,-1.2566907276473498),
     "island_equal":(5.383529291134088,-1.3647844594338443)},
  "measured_3p2r2s":{
     "cell_equal":(1.0160798695305384,-0.1481462709020268),
     "island_equal":(2.938700414161954,1.8406806075457232)}
}

def source_components(source_island_cells,original_species,taxon_ref):
    """Same original photographed source cells; no species/colour imputation."""
    source=original_species.loc[original_species.target_id.isin(
        set(source_island_cells.target_id))].copy()
    panel.require(len(source)>0 and not source.duplicated([
        "target_id",taxa.SPECIES_COL]).any(),
        "Source taxon-photo-cell rows are duplicated or missing")
    part=taxon_ref.loc[:,["target_id","taxon_average_white_fraction",
                         "source_singleton_taxa","taxon_species_count"]].copy()
    a=source_island_cells.merge(part,on="target_id",how="left",
                                validate="one_to_one")
    n=source.groupby("target_id")[taxa.SPECIES_COL].nunique()
    panel.require(a.taxon_average_white_fraction.notna().all() and
                  a.target_id.isin(n.index).all() and
                  np.array_equal(a.set_index("target_id").loc[n.index,
                      "accepted_species_n"].to_numpy(int),n.to_numpy(int)),
                  "Original actual WCVP photo taxa do not reproduce cell species")
    panel.require(a.taxon_species_count.eq(a.accepted_species_n).all(),
                  "Unverified species-equal photographic reference")
    a["taxon_reference_white_pp"]=100*a.taxon_average_white_fraction
    a["taxon_relative_photo_pp"]=a.white_pp-a.taxon_reference_white_pp
    panel.require(np.allclose(
        a.white_pp,a.taxon_reference_white_pp+a.taxon_relative_photo_pp,
        atol=1e-9,rtol=0),"Observed ROI outcome decomposition not exact")
    a["observed_taxa_in_cell"]=a.taxon_species_count
    return a

def fitted_term_differences(frame,col,area_centre,distance_centre,equal_island):
    """Fit all responses using mathematically identical X, weights, and islands."""
    d=frame.copy()
    d["log10_area_c"]=np.log10(d.area_km2)-area_centre
    d["log10_distance_c"]=np.log10(1+d.distance_mainland_vertex_km)-distance_centre
    panel.require(d.land_id.nunique()>=10 and d.year.nunique()>=2 and
                  d[["white_pp","log_light","log10_area_c","log10_distance_c",
                     "bio5","log_rain","log_pop","elevation","log_attempts",
                     col]].notna().all().all(),
                  "Missing original source outcome/controls or insufficient islands")
    formula=(col+" ~ log_light + log_light:log10_area_c + "
             "log_light:log10_distance_c + "+BASE_COVARS+
             " + C(land_id) + C(year)")
    if equal_island:
        weights=1/d.groupby("land_id",observed=True).land_id.transform(
             "size").astype(float)
        panel.require(np.allclose(weights.groupby(d.land_id).sum(),1),
                      "Islands do not have equal total source photo weight")
        fit=smf.wls(formula,d,weights=weights).fit(
            cov_type="cluster",cov_kwds={"groups":d.land_id},use_t=True)
    else:
        fit=smf.ols(formula,d).fit(
            cov_type="cluster",cov_kwds={"groups":d.land_id},use_t=True)
    terms={}
    for term in TERMS:
        panel.require(term in fit.params,"Source original term not identified")
        lo,hi=fit.conf_int().loc[term]
        nums={"beta":float(fit.params[term]),"SE":float(fit.bse[term]),
              "ci_low":float(lo),"ci_high":float(hi),
              "p":float(fit.pvalues[term])}
        panel.require(np.isfinite(list(nums.values())).all(),
                      "Island area-distance fit had nonfinite terms")
        terms[term]=nums
    return {"terms":terms,"photo_cell_years":len(d),
            "unique_islands":int(d.land_id.nunique()),
            "formula":formula,
            "weighting":"equal_island" if equal_island else "equal_cell_year"}

def exact_decomposition_model(data,observed,reference,deviation,
                              area_centre,distance_centre,equal_island):
    """Coefficient identity on the EXACT same rows (not causal contribution)."""
    y=data[observed].to_numpy(float)
    r=data[reference].to_numpy(float)
    v=data[deviation].to_numpy(float)
    panel.require(np.allclose(y,r+v,atol=1e-9,rtol=0),
                  "Attempt to decompose different photographed endpoints")
    outcomes={}
    for label,col in (("observed",observed),("taxon_photo_reference",reference),
                      ("taxon_relative_photo_deviation",deviation)):
        outcomes[label]=fitted_term_differences(
            data,col,area_centre,distance_centre,equal_island)
    for term in TERMS:
        observed_beta=outcomes["observed"]["terms"][term]["beta"]
        split_sum=(outcomes["taxon_photo_reference"]["terms"][term]["beta"]+
                   outcomes["taxon_relative_photo_deviation"]["terms"][term]["beta"])
        panel.require(abs(observed_beta-split_sum)<1e-6,
                      "Coefficient split is not exact for "+term)
    return {"results":outcomes,"both_interactions_reconstruct_exactly":True,
            "within_source_reference_not_genetic_or_independent_trait":True,
            "components_are_descriptive_not_mediation_proportions":True}

def earlier_year_reference_holdout(species,full_island):
    """Photographically observed prior taxa only; do NOT invent absent taxa."""
    prior=species.loc[species.year.isin(PRIOR)]
    reference=prior.groupby(taxa.SPECIES_COL,observed=True).agg(
        historical_white=("photo_white_fraction","mean"),
        historical_observed_cells=("target_id","nunique")).reset_index()
    future=species.loc[species.year.isin(HOLDOUT)&
                       species.target_id.isin(set(full_island.target_id))]
    observed=future.merge(reference,on=taxa.SPECIES_COL,how="inner",
                           validate="many_to_one")
    panel.require(observed.year.isin(HOLDOUT).all() and
                  observed.historical_white.between(0,1).all() and
                  observed.photo_white_fraction.between(0,1).all(),
                  "Holdout future photography leaked into earlier-year reference")
    cell=observed.groupby("target_id",observed=True).agg(
         n_observed_reference_taxa=(taxa.SPECIES_COL,"nunique"),
         prior_photo_reference_pp=("historical_white",lambda a:100*a.mean()),
         current_photo_white_pp=("photo_white_fraction",lambda a:100*a.mean())
         ).reset_index()
    cell["heldout_taxon_photo_delta_pp"]=(
        cell.current_photo_white_pp-cell.prior_photo_reference_pp)
    cell=cell.loc[cell.n_observed_reference_taxa.ge(2)].copy()
    joined=full_island.loc[full_island.year.isin(HOLDOUT)].merge(
        cell,on="target_id",how="inner",validate="one_to_one")
    panel.require(joined.n_observed_reference_taxa.le(
        joined.accepted_species_n).all(),
        "Holdout reference includes unobserved accepted species")
    panel.require(np.allclose(
        joined.current_photo_white_pp,
        joined.prior_photo_reference_pp+joined.heldout_taxon_photo_delta_pp,
        atol=1e-8,rtol=0),
        "Heldout current species source ROI not exactly decomposable")
    # Need actual different geographic sites within each island, not merely
    # two years for the same photographic 0.05° grid cell.
    source=joined.target_id.astype(str).str.split("|",regex=False,expand=True)
    panel.require(source.shape[1]==3 and
                  source[2].eq(joined.land_id.astype(str)).all(),
                  "Holdout GSHHG landmass identity drift")
    joined["site_id"]=source[1].to_numpy()
    by_island=joined.groupby("land_id").agg(
       real_different_sites=("site_id","nunique"),
       genuine_exposure_values=("log_light","nunique"))
    ids=by_island.index[(by_island.real_different_sites>=2)&
                        (by_island.genuine_exposure_values>=2)]
    supported=joined.loc[joined.land_id.isin(ids)].copy()
    return supported,{
       "prior_source_years":list(PRIOR),
       "untouched_test_observation_years":list(HOLDOUT),
       "source_reference_taxa_with_prior_photos":len(reference),
       "test_cell_years_having_2plus_genuinely_prior_photographed_species":len(joined),
       "test_islands_having_2plus_true_different_sites":int(len(ids)),
       "test_selected_cell_years":len(supported),
       "n_observed_reference_taxon_cell_years":int(
           supported.n_observed_reference_taxa.sum()) if len(supported) else 0,
       "holdout_species_colours_not_necessarily_same_genotypes":True,
       "holdout_photo_years_not_used_to_define_historical_taxon_reference":True,
       "no_missing_species_or_colour_imputation":True
    }

def run(paths,coast,out):
    out=Path(out)
    panel.require(not out.exists(),"Do not overwrite ALAN source evidence")
    frames=tier.build_frozen_five_years(paths)
    cohorts=tier.by_year_tier(frames)
    canonical,wcvp=panel.load_frozen_photos(paths)
    broad=cohorts["measured_3p2r2s"]
    panel.require(len(broad)==20543,"Original true photo-source broad cohort drift")
    species=taxa.real_species_cell_years(paths,canonical,wcvp,broad)
    source_taxa,taxonref=taxa.source_taxon_reference_white_components(species)
    land=isolation.load_signed_land(coast)
    meta=land[["land_id","land_class","area_km2"]]
    selected={}
    for name,n in COHORT_N.items():
        frame=cohorts[name]
        panel.require(len(frame)==n and frame.target_id.is_unique and
                      set(frame.year)==set(YEARS),
                      "Source photo actual cohort/year drift")
        all_island,within,groups=area.join_area_and_select(frame,meta)
        panel.require(len(within)==SOURCE_CELL_N[name] and
                      within.land_id.nunique()==SOURCE_ISLAND_N[name],
                      "Previously signed multi-site island sample altered")
        selected[name]=(all_island,within,groups)
    source_ids=set().union(*(set(w.land_id) for _,w,_ in selected.values()))
    distances=isolation.isolated_mainland_distance_km(land,source_ids)
    core=isolation.complete_island_frame(selected["strict_5p3r3s"][1],distances)
    center_area=float(np.log10(core.groupby("land_id").area_km2.first()).median())
    center_dist=float(np.log10(1+core.groupby("land_id")
                 .distance_mainland_vertex_km.first()).median())
    results={}
    for name,(all_island,within,_) in selected.items():
        photo=isolation.complete_island_frame(within,distances)
        analysis=source_components(photo,species,taxonref)
        report={
          "real_photo_cell_years":len(analysis),
          "real_islands":int(analysis.land_id.nunique()),
          "real_accepted_taxon_cell_years":int(analysis.observed_taxa_in_cell.sum()),
          "same_archive_singleton_taxon_cell_years":int(
             analysis.source_singleton_taxa.sum()),
          "same_archive_taxon_reference_training_cells":int(len(broad)),
          "decomposition":{},
        }
        for weight,eq in (("cell_equal",False),("island_equal",True)):
            d=exact_decomposition_model(
                analysis,"white_pp","taxon_reference_white_pp",
                "taxon_relative_photo_pp",center_area,center_dist,eq)
            old=isolation.ORIGINAL_BETAS if hasattr(isolation,"ORIGINAL_BETAS") else None
            # Use pinned, previously published ALAN joint numeric results.
            expected={
              "strict_5p3r3s":{
                "cell_equal":(2.216703936269394,-1.2566907276473498),
                "island_equal":(5.383529291134088,-1.3647844594338443)},
              "measured_3p2r2s":{
                "cell_equal":(1.0160798695305384,-0.1481462709020268),
                "island_equal":(2.938700414161954,1.8406806075457232)}}[name][weight]
            for term,v in zip(TERMS,expected):
                got=d["results"]["observed"]["terms"][term]["beta"]
                panel.require(abs(got-v)<1e-5,
                     "Original island area/distance real photo regression was changed")
            report["decomposition"][weight]=d
        heldout,holdout_support=earlier_year_reference_holdout(species,photo)
        report["heldout_source_support"]=holdout_support
        # Holdout has a separate genuine shared-taxon target population.
        # Avoid a numerical model when different-site evidence is very sparse.
        eligible=(len(heldout)>=100 and heldout.land_id.nunique()>=20
                  and heldout.year.nunique()==2
                  and heldout.groupby("land_id").site_id.nunique().ge(2).all())
        report["heldout_source_support"]["estimable"]=bool(eligible)
        if eligible:
            report["heldout_different_population_models"]={}
            for weight,eq in (("cell_equal",False),("island_equal",True)):
                report["heldout_different_population_models"][weight]=(
                    exact_decomposition_model(
                      heldout,"current_photo_white_pp",
                      "prior_photo_reference_pp","heldout_taxon_photo_delta_pp",
                      center_area,center_dist,eq))
        else:
            report["heldout_source_support"]["reason"]=(
                 "insufficient genuinely repeated original test-year island "
                 "sites with >=2 independently previously photographed taxa")
        results[name]=report
    summary={
      "schema":"ALAN_real_island_area_nightlight_photo_taxon_decomposition_v1",
      "years":list(YEARS),
      "original_strict_photo_cell_years":7982,
      "original_broad_photo_cell_years":20543,
      "species_reference_source":"same original WCVP-verified photo archive mean across broad cohort; not independently measured stable colour",
      "accepted_species_actual_source_count":int(len(source_taxa)),
      "actual_accepted_species_photo_cell_years":int(len(species)),
      "area_definition":"signed GSHHG original EPSG:6933 polygon area, per tenfold area",
      "mainland_distance_definition":"signed GSHHG WGS84 outer shoreline vertex geodesic approximation; per tenfold (1+km)",
      "fit_model":"photo_white_pp ~ log_light + log_light:log10_area_c + log_light:log10_distance_c + log_pop + bio5 + log_rain + elevation + log_attempts + island FE + year FE",
      "comparison_estimand":"same strict/broad original 44/95 islands and 650/1884 photo cell-years",
      "target_year_holdout_ref_uses_no_target_year_photographs":True,
      "island_equal_weights_one_per_GSHHG_island":True,
      "taxon_reference_components_are_algebraic_not_causal_turnover":True,
      "original_images_downloaded_or_reinferred":0,
      "no_missing_species_or_hue_imputation":True,
      "results":results,
      "source_SHA256":{**{k:panel.sha(v) for k,v in paths.items()},
                       "GSHHG_original_archive":isolation.COAST_SHA}
    }
    out.mkdir(parents=True)
    (out/"ALAN_island_area_accepted_taxon_photo_decomposition_aggregate.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("ALAN_REAL_ISLAND_AREA_TAXON_COMPONENTS",json.dumps({
        k:{
         "islands":v["real_islands"],
         "cell_years":v["real_photo_cell_years"],
         "source_accepted_taxon_cell_years":v["real_accepted_taxon_cell_years"],
         "heldout_support":v["heldout_source_support"],
         "models":{
           w:{component:{
               term:{kk:z[kk] for kk in ("beta","ci_low","ci_high","p")}
               for term,z in v["decomposition"][w]["results"][component]["terms"].items()}
               for component in ("observed","taxon_photo_reference",
                                 "taxon_relative_photo_deviation")}
             for w in ("cell_equal","island_equal")}
        } for k,v in results.items()},sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    keys=("canonical","wcvp",
          "geo2019","viirs2019","pop2019",
          "geo2021","viirs2021","pop2021",
          "geo2023","viirs2023","pop2023",
          "geo2024","cn2024","uc2024","uc2024receipt",
          "geo2025","geo2025receipt","viirs2025")
    for key in keys:
        parser.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    parser.add_argument("--coast",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    run({k:getattr(a,k) for k in keys},a.coast,a.out)
