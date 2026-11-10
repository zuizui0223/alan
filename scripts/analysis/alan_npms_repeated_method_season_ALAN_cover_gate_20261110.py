#!/usr/bin/env python3
"""NPMS original 2015–25 vegetation + genuine prior PCNL: source eligibility audit.

Use SAME real NPMS plot ID, true survey ID, seasonal window and source cover
scheme. Exact source-positive plant records only, no simulated zero occupancy.
GloPL is an independent species-level supplementation experiment, NOT PL in
an NPMS plot. Original 2014–24 calibrated PCNL is strictly previous year.
This script DOES NOT fit a biological response and does not select results by
p value. Logs only aggregate support, no volunteer ID or precise coordinates.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_npms_glopl_actual_species_repeat_support_20261110 import (
    data_read,names,exact_binomial,NPMS_MEMBERS)
from scripts.analysis.alan_npms_11yr_true_plot_priorPCNL_source_20261110 import (
    SCHEMA as LIGHT_SCHEMA)

SCHEMA="ALAN_actual_NPMS_sameplot_survey_season_cover_priorPCNL_support_v1"
EXPECTED_LIGHT_SHA="5ecae6a9b62de636c8e68fd694c69698eec154e5a8a0ba641beb5a7e6061e427"
DARK_MIN=.02

def require(ok,message):
    if not bool(ok):raise ValueError("ALAN_NPMS_POSITIVE_COVER_SOURCE: "+str(message))

def checksum(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for q in iter(lambda:f.read(4<<20),b""):h.update(q)
    return h.hexdigest()

def true_raw_cover_codes(z):
    key="data/dominscores.csv"
    with zipfile.ZipFile(z) as a:
        source=a.getinfo(key)
        size,crc=NPMS_MEMBERS[key]
        require(source.file_size==size and format(source.CRC,"08x")==crc,
                "Published DOMIN scores source has changed")
        with a.open(key) as f:
            d=pd.read_csv(f,dtype=str,keep_default_na=False)
    require({"domin","dominIndex","midPoint"}.issubset(d.columns),
            "No original cover source code mapping")
    return d

def source_frame(npms,glopl,light,receipt):
    r=json.loads(Path(receipt).read_text())
    require(r["schema"]==LIGHT_SCHEMA and
            r["genuine_source_2015_2025_original_plot_years"]==16960 and
            r["genuine_source_distinct_original_NPMS_plots"]==5360 and
            r["original_plotyear_radiance_source_CSV_SHA256"]==EXPECTED_LIGHT_SHA and
            checksum(light)==EXPECTED_LIGHT_SHA,
            "Original NPMS 16,960 site-year historical PCNL source hash mismatch")
    src=pd.read_csv(light,dtype={"location_id":str},low_memory=False)
    require(len(src)==16960 and not src.duplicated(["location_id","year"]).any() and
            (src.prior_PCNL_year==src.year-1).all(),
            "Original last-year PCNL site/year join changed")
    tbl,pl=data_read(npms,glopl)
    occ,samp=tbl["occurrences"],tbl["sampleinfo"]
    loc=tbl["locations"]
    require(samp.id.is_unique and loc.location_id.is_unique and
            {"survey_id","location_id","id","date_start"}.issubset(samp.columns),
            "Unknown original plot and survey ID")
    samp=samp.loc[samp.location_id.isin(set(loc.location_id))].copy()
    samp["date"]=pd.to_datetime(samp.date_start,utc=True,errors="coerce")
    samp["year"]=samp.date.dt.year
    samp["month"]=samp.date.dt.month
    samp["doy"]=samp.date.dt.dayofyear
    samp["season"]=np.select([samp.month.isin([4,5,6]),samp.month.isin([7,8,9])],
                             ["spring","summer"],default="other")
    samp["survey_method_id"]=samp.survey_id.astype(str).str.strip()
    samp["source_plot_type"]=samp.location_id.map(
        loc.set_index("location_id").location_type)
    samp=samp.merge(src[["location_id","year","prior_PCNL_year",
        "source_true_prior_PCNl","source_prior_log1p_PCNL",
        "site_is_spatially_fixed_within_500metres"]],
        on=["location_id","year"],how="left",validate="many_to_one",indicator=True)
    require(samp._merge.eq("both").all() and (samp.prior_PCNL_year==samp.year-1).all(),
            "Original measured NPMS sample and annual raster source keys altered")
    samp=samp.drop(columns="_merge")
    cover=true_raw_cover_codes(npms)
    o=occ[["sample_id","preferred_taxon","record_status","domin"]].copy()
    o["species"]=names(o.preferred_taxon)
    o=o.loc[exact_binomial(o.species)].copy()
    o=o.merge(samp[["id","location_id","year","doy","month","season","survey_method_id",
                    "source_plot_type","source_prior_log1p_PCNL","prior_PCNL_year"]],
             left_on="sample_id",right_on="id",how="inner",validate="many_to_one")
    o["source_original_domin"]=o.domin.astype(str).str.strip()
    o["original_record_status"]=o.record_status.astype(str).str.strip()
    # Preserve genuinely distinct source statuses and schema encodings.
    o["domin_encoding"]=np.select([
       o.source_original_domin.str.fullmatch(r"\d+"),
       o.source_original_domin.str.contains("%",regex=False),
       o.source_original_domin.eq("")],
       ["numeric_source_code","percent_text_range","missing"],
       default="other_source_text")
    ref=cover[["domin","dominIndex","midPoint"]].copy()
    ref["domin"]=ref.domin.astype(str).str.strip()
    ref["cover_midpoint_from_original_source"]=pd.to_numeric(
        ref.midPoint,errors="coerce")
    # A cover code is considered source-mappable ONLY if original source
    # dictionary has one consistent finite cover midpoint in [0,100].
    original_mappings=ref.groupby("domin",sort=False).agg(
        n_distinct_midpoints=("cover_midpoint_from_original_source","nunique"),
        midpoint=("cover_midpoint_from_original_source","first"),
        lookup_rows=("domin","size")).reset_index()
    valid_mappings=original_mappings.loc[
        original_mappings.n_distinct_midpoints.eq(1)&
        original_mappings.midpoint.between(0,1)].copy()
    require(not valid_mappings.domin.duplicated().any(),
            "Original DOMIN coding lookup not unique")
    o=o.merge(valid_mappings[["domin","midpoint"]],left_on="source_original_domin",
              right_on="domin",how="left",validate="many_to_one",suffixes=("","_lookup"))
    o["cover_midpoint_from_original_source"]=o.midpoint
    # Report SOURCE provenance by scheme; do NOT pool numeric source values
    # from incompatible survey methods unless exact original code lookup
    # demonstrates identical measured percent unit.
    o["log1p_cover"]=np.log1p(o.cover_midpoint_from_original_source)
    pl["species"]=names(pl.Species_accepted_names.str.replace("_"," ",regex=False))
    pl=pl.loc[exact_binomial(pl.species)].copy()
    pl["study_key"]=np.where(pl.DOI.astype(str).str.strip().ne(""),
        names(pl.DOI).str.casefold(),names(pl.Author)+"|"+names(pl.Year))
    glopl_studies=pl.groupby("species").study_key.nunique()
    o["has_original_GloPL_case"]=o.species.isin(set(glopl_studies.index))
    o["has_two_independent_original_GloPL_studies"]=o.species.isin(
        set(glopl_studies.index[glopl_studies.ge(2)]))
    return o,samp,cover,original_mappings,src,r

def run(npms,glopl,light,receipt,out):
    o,samp,cover,mapping,source,rec=source_frame(npms,glopl,light,receipt)
    same=samp.loc[samp.source_prior_log1p_PCNL.notna()&
                  samp.season.isin(["spring","summer"])].copy()
    by=("location_id","survey_method_id","season")
    groups=same.groupby(list(by),sort=False).agg(
        years=("year","nunique"),
        original_visits=("id","nunique"),
        min_light=("source_prior_log1p_PCNL","min"),
        max_light=("source_prior_log1p_PCNL","max"),
        n_light_values=("source_prior_log1p_PCNL","nunique"))
    allowed=groups.loc[groups.years.ge(2)&
                       (groups.max_light-groups.min_light>=DARK_MIN)&
                       groups.n_light_values.ge(2)].copy()
    eligible=allowed.index
    source_positive=o.loc[o.source_prior_log1p_PCNL.notna()&
                         o.season.isin(["spring","summer"])].copy()
    specieskey=list(by)+["species"]
    positive=source_positive.drop_duplicates(specieskey+["year","sample_id"])
    source_support=positive.groupby(specieskey,sort=False).agg(
       source_positive_years=("year","nunique"),
       source_visited_rows=("sample_id","nunique"),
       n_numeric_cover=("cover_midpoint_from_original_source",lambda s:int(s.notna().sum())),
       measured_cover_years=("year",lambda s:int(s.nunique())))
    positive_repeated=source_support.loc[source_support.source_positive_years.ge(2)]
    source_positive=source_positive.merge(allowed.reset_index()[list(by)],
                         on=list(by),how="inner",validate="many_to_one")
    measured=source_positive.loc[
        source_positive.cover_midpoint_from_original_source.notna()].copy()
    agg=measured.groupby(specieskey+["year"],sort=False).agg(
        median_cover=("cover_midpoint_from_original_source","median"),
        median_source_day_of_year=("doy","median"),
        n_visits_with_cover=("sample_id","nunique"),
        source_domin_encodings=("domin_encoding",lambda x:";".join(sorted(set(x)))),
        source_record_status=("original_record_status",lambda x:";".join(sorted(set(x)))),
        original_annual_priorPCNL=("source_prior_log1p_PCNL","first"),
        source_prior_unique_values=("source_prior_log1p_PCNL","nunique")).reset_index()
    require(agg.source_prior_unique_values.le(1).all(),
            "Within plot×year original year-prior source light must be identical")
    direct=agg.groupby(specieskey,sort=False).agg(
        years=("year","nunique"),
        n_numeric_cover_rows=("median_cover","size"),
        n_pcNl_values=("original_annual_priorPCNL","nunique"),
        max_light=("original_annual_priorPCNL","max"),
        min_light=("original_annual_priorPCNL","min"),
        n_cover_source_encodings=("source_domin_encodings","nunique")).reset_index()
    direct=direct.loc[direct.years.ge(2)&
                      (direct.max_light-direct.min_light>=DARK_MIN)&
                      direct.n_pcNl_values.ge(2)]
    multi=set(source_positive.loc[
        source_positive.has_two_independent_original_GloPL_studies,"species"])
    summary={
      "schema":"ALAN_actual_NPMS_sameplot_survey_season_cover_priorPCNL_support_v1",
      "source_NPMS_EIDC_DOI":"10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8",
      "source_prior_PCNL_DOI":"10.5281/zenodo.17013414",
      "original_16960_plotyear_light_SHA256":EXPECTED_LIGHT_SHA,
      "source_NPMS_distinct_plots":int(samp.location_id.nunique()),
      "original_NPMS_valid_actual_sample_visits":len(samp),
      "source_2015to2025_same_surveyID_season_plot_groups_with_two_years_and_prior_light_span_ge0p02":len(allowed),
      "source_2015to2025_same_surveyID_season_eligible_unique_plots":len(set(i[0] for i in eligible)),
      "source_same_method_season_prior_PCNL_logspan_quantiles":{
         str(q):float((allowed.max_light-allowed.min_light).quantile(q))
         for q in (.1,.25,.50,.75,.9)} if len(allowed) else {},
      "source_positive_NPMS_named_plant_records_all_source_methods":len(o),
      "source_positive_plant_records_with_numeric_cover_lookup":int(
          o.cover_midpoint_from_original_source.notna().sum()),
      "source_positive_plant_domin_code_types_by_original_surveyID":{
        str(s):{str(k):int(v) for k,v in x.domin_encoding.value_counts().items()}
        for s,x in o.groupby("survey_method_id")},
      "source_positive_plant_cover_lookup_resolved_by_original_surveyID":{
        str(s):{"n_records":len(x),
                "cover_midpoint_resolution_fraction":float(
                   x.cover_midpoint_from_original_source.notna().mean()),
                "domin_examples":{str(k):int(v) for k,v in
                   x.source_original_domin.value_counts().head(16).items()}}
        for s,x in o.groupby("survey_method_id")},
      "source_CRC_pinned_DOMIN_midPoint_unit":"fraction of cover [0,1], not percent [0,100]",
      "source_records_have_actual_doy_for_within_season_phenology_adjustment":True,
      "source_original_DOMIN_score_table":cover[["domin","dominIndex","midPoint"]].to_dict("records"),
      "source_original_DOMIN_unique_code_midpoints":
          mapping[["domin","midpoint","n_distinct_midpoints","lookup_rows"]].to_dict("records"),
      "source_GloPL_species_with_exact_name_in_NPMS":79,
      "source_valid_NPMS_same_method_season_light_changing_plant_positive_rows":
         len(source_positive),
      "source_valid_NPMS_same_method_season_light_changing_plant_positive_species":
         int(source_positive.species.nunique()),
      "source_valid_NPMS_same_method_season_light_changing_plant_positive_plots":
         int(source_positive.location_id.nunique()),
      "source_valid_NPMS_same_method_season_light_changing_plant_positive_GloPL_taxa":
         int(source_positive.loc[source_positive.has_original_GloPL_case,"species"].nunique()),
      "source_valid_NPMS_same_method_season_light_changing_plant_positive_multi_study_GloPL_taxa":
         int(source_positive.loc[source_positive.has_two_independent_original_GloPL_studies,"species"].nunique()),
      "source_positive_NPMS_cover_plot_species_method_season_repeated_and_light_change":
         len(direct),
      "source_positive_NPMS_cover_eligible_distinct_plots":
         int(direct.location_id.nunique()),
      "source_positive_NPMS_cover_eligible_distinct_species":
         int(direct.species.nunique()),
      "source_positive_NPMS_cover_eligible_distinct_PCNl_plot_year_rows":
         len(agg),
      "source_positive_NPMS_cover_eligible_multi_study_GloPL_plot_species_groups":
         int(direct.species.isin(multi).sum()),
      "source_not_recorded_plant_never_imputed_absent":True,
      "positive_only_cover_conditioned_on_detection":True,
      "source_DOMIN_numeric_code_mapping_requires_published_lookup_not_implicit_1to10":True,
      "no_causal_ALAN_plant_occupancy_or_effective_moth_pollination_claimed":True}
    out=Path(out);require(not out.exists(),"No overwrite source eligibility audit")
    out.mkdir(parents=True)
    # Export only eligible source-positive plot×taxon×method×season×year
    # numeric covers, no observer or precise plot coordinates.
    agg.to_csv(out/"NPMS_source_positive_comparable_plant_cover_plot_method_season_year.csv.gz",
               index=False,compression="gzip")
    (out/"NPMS_2015to2025_repeated_plot_method_season_alanchanging_cover_source_receipt.json").write_text(
        json.dumps(summary,sort_keys=True,indent=2,default=str)+"\n")
    print("ALAN_TRUE_NPMS_SAME_METHOD_SEASON_ALAN_COVER_SUPPORT",
          json.dumps(summary,sort_keys=True,default=str),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser()
    for name in ("npms","glopl","light","receipt","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.npms,a.glopl,a.light,a.receipt,a.out)

if __name__=="__main__":main()
