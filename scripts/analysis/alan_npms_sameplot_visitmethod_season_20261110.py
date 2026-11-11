#!/usr/bin/env python3
"""NPMS 2015-25 true repeat-plot, visit-method and season comparability audit.

This is a *source-eligibility* measurement for future time-varying ALAN vs
documented vegetation change, NOT a plant disappearance, occupancy,
effective pollination or pollen limitation estimate.

Source: published NPMS EIDC 10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8,
seven original CSV size+CRC gates; published GloPL exact SHA.
Never fill nonrecorded species with biological zero. Preserve same actual
plot ID, sample/year, survey IDs and season; empirical comparability
gates must be measured before claiming within-plot change.
"""
from __future__ import annotations

import argparse
import collections
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.alan_npms_glopl_actual_species_repeat_support_20261110 import (
    data_read, names, exact_binomial, GLOPL_SHA)

SCHEMA="ALAN_real_NPMS_2015_2025_plot_season_method_source_support_v1"

def require(ok,reason):
    if not bool(ok): raise ValueError("ALAN_NPMS_REPEAT_DESIGN: "+str(reason))

def as_counts(s,limit=25):
    return {str(k):int(v) for k,v in s.fillna("UNKNOWN").astype(str).value_counts().head(limit).items()}

def season_for_month(m):
    if m in (4,5,6):return "spring_AprJun"
    if m in (7,8,9):return "summer_JulSep"
    return "off_season_other"

def source_vegetation_design(package,glopl):
    t,g=data_read(package,glopl)
    sample=t["sampleinfo"].copy()
    occ=t["occurrences"].copy()
    loc=t["locations"].copy()
    attr=t["attributes"].copy()
    require({"id","survey_id","date_start","location_id",
             "LATITUDE","LONGITUDE","country","created_by_id"}<=set(sample),
            "Original NPMS survey-unit metadata changed")
    require({"sample_id","preferred_taxon","record_status","domin"}<=set(occ),
            "Original documented species table changed")
    require({"location_id","location_type"}<=set(loc) and
            {"sample_id","caption","text_value"}<=set(attr),
            "Original plot method or sampled habitat source missing")
    require(len(sample)==27611 and len(occ)==271701 and
            len(loc)==5360 and len(attr)==293561 and len(g)==2969,
            "Original NPMS or GloPL source row frame changed")
    require(sample.id.is_unique and loc.location_id.is_unique,
            "Repeated original sample ID or GSHHG-independent NPMS plot ID")
    sample["visit_date"]=pd.to_datetime(sample.date_start,errors="coerce",utc=True)
    sample["year"]=sample.visit_date.dt.year
    sample["month"]=sample.visit_date.dt.month
    sample["doy"]=sample.visit_date.dt.dayofyear
    sample["season"]=sample.month.map(season_for_month)
    sample["lat"]=pd.to_numeric(sample.LATITUDE,errors="coerce")
    sample["lon"]=pd.to_numeric(sample.LONGITUDE,errors="coerce")
    sample["location_id"]=sample.location_id.astype(str)
    valid_location=sample.location_id.isin(set(loc.location_id))
    valid_geo=sample.lat.between(-90,90)&sample.lon.between(-180,180)
    valid_date=sample.year.between(2015,2025)
    matched=sample.loc[valid_location&valid_geo&valid_date].copy()
    require(len(matched)==27603 and matched.location_id.nunique()==5360,
            "The previously signed NPMS georeferenced sample frame changed")
    matched=matched.merge(loc[["location_id","location_type"]],
                          on="location_id",how="left",validate="many_to_one")
    require(matched.location_type.notna().all(),"Genuine plot location method unresolved")
    require(not matched.duplicated(["id"]).any(),"Visits no longer have original ID")

    # The sample attributes contain habitat + source method evidence, but a
    # caption may be present more than once on a sample. Do not pick an
    # arbitrary value: stable only if the unique nonblank one exists.
    target_attr=("NPMS Habitat","NPMS Survey 1 ID","Data Source","Data Source Version")
    all_attr=attr.loc[attr.caption.isin(target_attr)].copy()
    def single_nonempty(z):
        x=sorted(set(v.strip() for v in z.astype(str) if v.strip()))
        return x[0] if len(x)==1 else None
    attr_group=all_attr.groupby(["sample_id","caption"],sort=False).text_value.agg(
        single_nonempty).unstack("caption")
    visits=matched.join(attr_group,on="id")
    require(len(visits)==len(matched),"Original attribute pivot changed NPMS visit count")
    visits["survey_token"]=visits.survey_id.fillna("").astype(str).str.strip()
    visits["observer_token"]=visits.created_by_id.fillna("").astype(str).str.strip()
    visits["plot_observation_type"]=visits.location_type.fillna("").astype(str)
    visits["season_method_key"]=visits["season"]+"|"+visits["plot_observation_type"]
    visits["plot_season_key"]=visits.location_id+"|"+visits.season
    visits["plot_year_key"]=visits.location_id+"|"+visits.year.astype(int).astype(str)
    visits["plot_year_season_key"]=visits.plot_year_key+"|"+visits.season

    # Each actual visit is a sampling opportunity. Output only denominator
    # support, not inferred absence for plants never recorded in a plot-year.
    names_npms=names(occ.preferred_taxon)
    x=occ.loc[exact_binomial(names_npms)].copy()
    x["species"]=names_npms.loc[x.index]
    x=x.merge(visits[["id","location_id","year","month","doy","season",
                       "survey_token","plot_observation_type"]],
              left_on="sample_id",right_on="id",how="inner",
              validate="many_to_one")
    require(x.sample_id.nunique()>0 and len(x)>50000,
            "Original accepted actual plant taxon row source missing")

    gl=g.copy()
    gl["species"]=names(gl.Species_accepted_names.str.replace("_"," ",regex=False))
    gl=gl.loc[exact_binomial(gl.species)].copy()
    gl["study_key"]=np.where(gl.DOI.astype(str).str.strip().ne(""),
        names(gl.DOI).str.casefold(),names(gl.Author)+"|"+names(gl.Year))
    count_study=gl.groupby("species").study_key.nunique()
    taxa=set(gl.species.unique())
    replicated_study_taxa=set(count_study.index[count_study.ge(2)])
    matched_sp=x.loc[x.species.isin(taxa)].copy()
    measured_multi_study=x.loc[x.species.isin(replicated_study_taxa)].copy()

    # Source-positive-only site/taxon/year coverage. A nonrecord is UNKNOWN,
    # not a local occupancy 0.
    documented=matched_sp.drop_duplicates(["location_id","year","species"])
    repeat_taxa=(
        documented.groupby(["location_id","species"],sort=False).year.nunique())
    rep_two=repeat_taxa.loc[repeat_taxa.ge(2)]
    multi_study_pos=measured_multi_study.drop_duplicates(
        ["location_id","year","species"])
    multirep=multi_study_pos.groupby(["location_id","species"]).year.nunique()

    plotyears=visits.groupby("location_id",sort=False).year.nunique()
    plotseasonyears=visits.groupby(["location_id","season"],sort=False).year.nunique()
    plotseasons=visits.groupby(["location_id","year","season"],sort=False).id.nunique()
    doys=visits.groupby(["location_id","season"],sort=False).doy.agg(
        lambda a: float(a.max()-a.min()) if a.notna().all() else np.nan)
    obsyears=visits.groupby(["location_id","observer_token"],sort=False).year.nunique()
    surveyyears=visits.groupby(["location_id","survey_token"],sort=False).year.nunique()
    true_same_survey_and_season_years=visits.groupby(
        ["location_id","survey_token","season"],sort=False).year.nunique()
    strict_same_source_plots=set(
        true_same_survey_and_season_years.loc[
           true_same_survey_and_season_years.ge(2)].reset_index().location_id)
    positive_same_survey_season=matched_sp.drop_duplicates(
        ["location_id","survey_token","season","year","species"])
    repeated_positive_same_survey_season=positive_same_survey_season.groupby(
        ["location_id","survey_token","season","species"],sort=False).year.nunique()
    repeated_multi_study_same_survey_season=positive_same_survey_season.loc[
        positive_same_survey_season.species.isin(replicated_study_taxa)].groupby(
        ["location_id","survey_token","season","species"],sort=False).year.nunique()

    # Is survey_id stable within plot across years? This is an *audit*, not
    # an unsupported guess that survey_id equals survey METHOD.
    sample_summary={
      "original_sample_visits":len(sample),
      "valid_true_location_year_coordinate_visits":len(visits),
      "actual_plot_locations":int(visits.location_id.nunique()),
      "actual_original_distinct_years":int(visits.year.nunique()),
      "actual_plot_year_units":int(visits[["location_id","year"]].drop_duplicates().shape[0]),
      "source_plots_revisited_2plus_years":int(plotyears.ge(2).sum()),
      "source_plots_revisited_3plus_years":int(plotyears.ge(3).sum()),
      "source_plot_season_groups_revisited_2plus_years":int(plotseasonyears.ge(2).sum()),
      "source_distinct_plots_revisited_in_same_season_2plus_years":int(
          plotseasonyears.loc[plotseasonyears.ge(2)].reset_index().location_id.nunique()),
      "source_plot_year_season_units":int(len(plotseasons)),
      "source_plot_year_season_units_with_multiple_visits":int(plotseasons.gt(1).sum()),
      "source_plot_season_revisit_doy_span_median":float(doys.dropna().median()),
      "source_plot_observer_tokens_revisited_2plus_years":int(obsyears.ge(2).sum()),
      "source_plots_with_any_repeat_survey_token_2plus_years":int(
          surveyyears.loc[surveyyears.ge(2)].reset_index().location_id.nunique()),
      "source_plots_with_same_survey_ID_AND_same_season_repeated_2plus_years":
          int(len(strict_same_source_plots)),
      "source_plot_survey_ID_season_groups_repeated_2plus_years":
          int(true_same_survey_and_season_years.ge(2).sum()),
      "original_survey_id_distinct_values":int(visits.survey_token.nunique()),
      "source_sample_survey_id_blank_visits":int(visits.survey_token.eq("").sum()),
      "source_unique_visitor_observer_ids":int(visits.observer_token.replace("",np.nan).nunique()),
      "source_plot_location_types":as_counts(visits.plot_observation_type),
      "source_sample_survey_id_top_categories":as_counts(visits.survey_token),
      "source_sample_season_counts":as_counts(visits.season),
      "source_sample_year_counts":{str(y):int(n) for y,n in visits.year.value_counts().sort_index().items()},
      "source_habitat_from_sampleattributes_present":int(visits["NPMS Habitat"].notna().sum()),
      "source_unique_sample_habitat_strings":int(visits["NPMS Habitat"].nunique()),
      "source_NPMS_habitat_top25":as_counts(visits["NPMS Habitat"],25),
      "source_data_source_top25":as_counts(visits["Data Source"],25),
      "source_plant_record_status_codes":as_counts(occ.record_status),
      "source_plant_domin_cover_codes":as_counts(occ.domin),
      "single_method_comparability_NOT_equivalent_to_similar_plot_type":True,
      "within_plot_month_alignment_NOT_asserted_from_season_only":True,
      "original_source_visit_no_plant_detection_cannot_automatically_mean_zero":True
    }
    plant_support={
      "original_unambiguous_binomial_plant_rows_linked_to_valid_visits":len(x),
      "genuine_GloPL_exact_source_species_in_NPMS":int(documented.species.nunique()),
      "original_NPMS_GloPL_site_year_species_recorded_positive":len(documented),
      "source_NPMS_plots_with_a_GloPL_species_recorded":int(documented.location_id.nunique()),
      "positive_recorded_plot_species_pairs_repeated_2plus_observation_years":int(len(rep_two)),
      "positive_repeat_2plus_GloPL_species":int(rep_two.reset_index().species.nunique()) if len(rep_two) else 0,
      "positive_repeat_2plus_GloPL_original_plots":int(rep_two.reset_index().location_id.nunique()) if len(rep_two) else 0,
      "original_GloPL_taxa_with_2plus_original_independent_studies_in_NPMS":
          int(multi_study_pos.species.nunique()),
      "original_multi_study_GloPL_positive_plot_year_species":len(multi_study_pos),
      "source_multi_study_species_plot_pairs_with_2plus_positive_years":int(multirep.ge(2).sum()),
      "source_same_survey_ID_AND_season_GloPL_positive_plot_species_repeat_groups":
           int(repeated_positive_same_survey_season.ge(2).sum()),
      "source_same_survey_ID_AND_season_independent_multi_study_GloPL_positive_repeat_groups":
           int(repeated_multi_study_same_survey_season.ge(2).sum()),
      "original_GloPL_species_PL_not_measured_at_NPMS_site":True,
      "source_photo_species_positive_repeat_does_not_establish_nonrecord_absence":True,
      "not_a_night_pollinator_or_moth_observation_dataset":True,
    }
    require(plant_support["genuine_GloPL_exact_source_species_in_NPMS"]==79 and
            plant_support["original_GloPL_taxa_with_2plus_original_independent_studies_in_NPMS"]==19,
            "Prior source evidence 79 and 19 species does not reconcile")
    return {"schema":SCHEMA,"NPMS_EIDC_original_DOI":
            "10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8",
            "GloPL_original_DOI":"10.5061/dryad.dt437",
            "GloPL_verified_source_SHA256":GLOPL_SHA,
            "original_NPMS_method_and_repeat_visit_support":sample_summary,
            "real_recorded_plant_species_with_independent_GloPL_support":plant_support,
            "effective_night_pollination_or_Moth_population_not_measured":True,
            "ALAN_emission_and_year_match_not_yet_fitted_for_NPMS_plots":True,
            "source_missing_species_zeros_not_imputed":True}

def run(package,glopl,out):
    out=Path(out)
    require(not out.exists(),"Refuse replace original repeated-plot audit")
    report=source_vegetation_design(package,glopl)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_REAL_NPMS_REVISIT_METHOD_SEASON_EFFORT_AUDIT",
          json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--package",type=Path,required=True)
    p.add_argument("--glopl",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.package,a.glopl,a.out)
