#!/usr/bin/env python3
"""Count genuine NPMS 2015–2025 repeated vegetation survey overlap with GloPL.

Names are matched as exact, unambiguous, full binomials only. Absence is NOT
inferred from missing NPMS species records, because NPMS field protocols may
be indicator lists rather than full vascular-plant censuses. GloPL
experimental pollen limitation is independent taxon-level background evidence,
never a measurement in a specific NPMS plot. All input sources are SHA-pinned.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile

import numpy as np
import pandas as pd

# EIDC dynamically generates ZIP wrappers: ZIP SHA may differ on re-request.
# Pin exact uncompressed CSV source member byte size and ZIP CRC from first audited
# 2015-2025 package; ZIP wrapper SHA is reported, never used as immutability proof.
NPMS_MEMBERS={
    "data/occurrences_2015to2025.csv":(18724281,"616d29d2"),
    "data/npmshabitatlookup_fine2broad.csv":(1411,"60a53016"),
    "data/sampleattributes_2015to2025.csv":(12159752,"9c9465a3"),
    "data/npmshabitatlookup_surv2broad.csv":(1848,"4b2a1f88"),
    "data/dominscores.csv":(433,"64c7fb11"),
    "data/locations_2015to2025.csv":(200780,"d5d18b79"),
    "data/sampleinfowithlatlong_2015to2025.csv":(4095902,"5e52f384"),
}
GLOPL_SHA="264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de"
KEYS={
    "occurrences":"data/occurrences_2015to2025.csv",
    "sampleinfo":"data/sampleinfowithlatlong_2015to2025.csv",
    "locations":"data/locations_2015to2025.csv",
    "attributes":"data/sampleattributes_2015to2025.csv",
}
SCHEMA="ALAN_true_NPMS_GloPL_repeat_plot_species_overlap_20261110_v1"
def require(ok,msg):
    if not ok:raise ValueError("ALAN_NPMS_GLOPL_OVERLAP: "+str(msg))
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()
def names(s):
    return s.fillna("").astype(str).str.strip().str.replace(r"\s+"," ",regex=True)
def exact_binomial(s):
    return names(s).str.fullmatch(r"[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+")
def data_read(package,glopl):
    require(sha(glopl)==GLOPL_SHA,"GloPL verified original Dryad source byte drift")
    out={}
    with zipfile.ZipFile(package) as z:
        require(z.testzip() is None,"Original NPMS EIDC ZIP failed CRC")
        entries={a.filename:a for a in z.infolist()}
        for name,(expected_bytes,crc) in NPMS_MEMBERS.items():
            require(name in entries and
                    entries[name].file_size==expected_bytes and
                    format(entries[name].CRC,"08x")==crc,
                    "NPMS original 2015-2025 uncompressed source member changed: "+name)
        require(set(KEYS.values()).issubset(z.namelist()),
                "Expected 2015–2025 original NPMS four source tables absent")
        for key,path in KEYS.items():
            with z.open(path) as stream:
                out[key]=pd.read_csv(stream,dtype=str,keep_default_na=False,
                                      encoding="utf-8-sig",low_memory=False)
    gl=pd.read_csv(glopl,dtype=str,keep_default_na=False,encoding="latin-1")
    require(len(out["occurrences"])==271701 and len(out["sampleinfo"])==27611 and
            len(out["locations"])==5360 and len(out["attributes"])==293561 and len(gl)==2969,
            "Original validated NPMS/GloPL source row totals changed")
    return out,gl
def run(package,glopl,out):
    tbl,gl=data_read(package,glopl)
    o=tbl["occurrences"].copy()
    sample=tbl["sampleinfo"].copy()
    loc=tbl["locations"].copy()
    attrs=tbl["attributes"].copy()
    require({"sample_id","preferred_taxon","record_status","domin"}<=set(o) and
            {"id","date_start","location_id","LATITUDE","LONGITUDE","country"}<=set(sample),
            "Original NPMS taxon, plot or date columns changed")
    require(sample.id.is_unique,"One NPMS sample identity repeated in original table")
    require(loc.location_id.is_unique,"Original NPMS location identifiers duplicated")
    sample["year"]=pd.to_datetime(sample.date_start,errors="coerce").dt.year
    sample["lat"]=pd.to_numeric(sample.LATITUDE,errors="coerce")
    sample["lon"]=pd.to_numeric(sample.LONGITUDE,errors="coerce")
    sample["site_id"]=sample.location_id
    valid_geo=sample.lat.between(-90,90)&sample.lon.between(-180,180)
    valid_year=sample.year.between(2015,2025)
    joined=o.merge(sample[["id","year","site_id","lat","lon","country"]],
                   left_on="sample_id",right_on="id",how="left",
                   validate="many_to_one",indicator=True)
    matched=joined._merge.eq("both")
    require(matched.mean()>.95,
            "NPMS occurrence-to-survey sample foreign-key linkage incomplete")
    # The original sample table includes visits whose location_id is blank.
    # These are NOT a reproducibly revisitable plot. Keep them in source audit,
    # but exclude them from plot×year×taxon coverage and temporal responses.
    valid_location_ids=set(loc.location_id.astype(str))
    # The NPMS extract uses the literal string 'NA' for an unassigned
    # location_id. A coordinate attached to a sample does not make 'NA'
    # a stable, independently revisitable plot. Require location foreign key.
    no_location=~sample.site_id.astype(str).isin(valid_location_ids)
    x=joined.loc[matched & joined.year.between(2015,2025) &
                 joined.site_id.astype(str).isin(valid_location_ids)].copy()
    x["accepted_name_candidate"]=names(x.preferred_taxon)
    binomial=exact_binomial(x.accepted_name_candidate)
    x=x.loc[binomial].copy()
    # Need true NPMS visit for documented recorded-species opportunities.
    site_species=x.groupby(["site_id","year","accepted_name_candidate"],observed=True).agg(
        n_source_occurrence_rows=("sample_id","size"),
        n_source_sample_visits=("sample_id","nunique")).reset_index()
    require(not site_species.duplicated(["site_id","year","accepted_name_candidate"]).any(),
            "NPMS source accepted species counted more than once per plot-year")
    original_GloPL_names=names(gl.Species_accepted_names)
    original_GloPL_name_format={
      "n_unique_original":int(original_GloPL_names.nunique()),
      "n_with_underscore":int(original_GloPL_names.str.contains("_",regex=False).sum()),
      "examples":original_GloPL_names.drop_duplicates().head(10).tolist()}
    # GloPL's original phylogenetic taxa names may be joined with a literal
    # underscore; replace it only at token separator, not by fuzzy synonym.
    gl["accepted_name_candidate"]=names(
        gl.Species_accepted_names.str.replace("_"," ",regex=False))
    gl=gl.loc[exact_binomial(gl.accepted_name_candidate)].copy()
    gl["PL_log_response_ratio"]=pd.to_numeric(gl.PL_Effect_Size,errors="coerce")
    gl["study_key"]=np.where(gl.DOI.astype(str).str.strip().ne(""),
        names(gl.DOI).str.casefold(),
        names(gl.Author)+"|"+names(gl.Year))
    require(gl.study_key.astype(str).str.strip().ne("").all(),
            "No independent experiment study key for pollen limitation")
    pl=gl.groupby("accepted_name_candidate",observed=True).agg(
        n_GloPL_cases=("study_key","size"),
        n_finite_PL_cases=("PL_log_response_ratio",lambda a: int(np.isfinite(a).sum())),
        n_GloPL_study_keys=("study_key","nunique"),
        median_GloPL_experiment_log_ratio=("PL_log_response_ratio","median")).reset_index()
    sp=site_species.groupby("accepted_name_candidate",observed=True).agg(
        n_NPMS_site_years=("year","size"),n_NPMS_distinct_sites=("site_id","nunique"),
        n_NPMS_recorded_occurrence_rows=("n_source_occurrence_rows","sum"),
        n_NPMS_recorded_survey_visits=("n_source_sample_visits","sum")).reset_index()
    taxon=sp.merge(pl,on="accepted_name_candidate",how="left",
                   validate="one_to_one",indicator=True)
    taxon["has_GloPL_experimental_supplementation_record"]=taxon._merge.eq("both")
    taxon=taxon.drop(columns="_merge")
    matched_names=set(taxon.loc[taxon.has_GloPL_experimental_supplementation_record,
                                 "accepted_name_candidate"])
    eligible=site_species.loc[site_species.accepted_name_candidate.isin(matched_names)]
    geo_sample=sample.loc[valid_geo & valid_year & ~no_location].copy()
    repeat=geo_sample.groupby("site_id",observed=True).year.nunique()
    location_types=loc.location_type.value_counts(dropna=False).to_dict()
    report={
      "schema":SCHEMA,
      "NPMS_source_zip_wrapper_SHA256_not_pinned":sha(package),
      "NPMS_pinned_original_CSV_members_byte_CRC":NPMS_MEMBERS, "GloPL_source_SHA256":GLOPL_SHA,
      "NPMS_2015_2025_source_plant_occurrence_rows":len(o),
      "NPMS_original_sample_visits":len(sample),
      "NPMS_original_location_rows":len(loc),
      "NPMS_source_sample_attributes":len(attrs),
      "NPMS_occurrences_with_matched_sample_keys":int(matched.sum()),
      "NPMS_source_sample_visits_with_valid_coordinate_and_year_and_site_id":len(geo_sample),
      "NPMS_sample_visits_with_missing_location_id":int(no_location.sum()),
      "NPMS_distinct_geocoded_sites":int(geo_sample.site_id.nunique()),
      "NPMS_sites_with_at_least_two_sample_years":int((repeat>=2).sum()),
      "NPMS_sites_with_at_least_three_sample_years":int((repeat>=3).sum()),
      "NPMS_binomial_source_species":len(sp),
      "NPMS_unique_site_year_species_units":len(site_species),
      "NPMS_source_sample_attribute_captions":attrs.caption.value_counts().head(25).to_dict(),
      "NPMS_source_record_status_counts":o.record_status.value_counts(dropna=False).to_dict(),
      "NPMS_source_location_types":location_types,
      "NPMS_source_country_sample_counts":sample.country.value_counts(dropna=False).to_dict(),
      "GloPL_original_unique_unambiguous_binomial_species":len(pl),
      "GloPL_original_name_format":original_GloPL_name_format,
      "NPMS_GloPL_exact_binomial_species_overlap":len(matched_names),
      "NPMS_site_year_species_units_with_any_GloPL_supplementation":len(eligible),
      "NPMS_unique_sites_with_GloPL_matched_species":int(eligible.site_id.nunique()),
      "NPMS_GloPL_matched_species_with_two_or_more_independent_studies":int(
          taxon.loc[taxon.has_GloPL_experimental_supplementation_record &
                     taxon.n_GloPL_study_keys.ge(2)].shape[0]),
      "island_GSHHG_unique_plot_assignments":None,
      "NPMS_repeat_plot_method_comparability_audited":False,
      "NPMS_non_detection_is_not_zero_presence":True,
      "GloPL_absent_species_PL_unknown_not_zero":True,
      "GloPL_species_median_PL_not_NPMS_plot_measured_PL":True,
      "GloPL_publication_Year_not_assumed_experiment_date":True,
      "validated_ALAN_light_by_NPMS_plant_plot":False,
      "pollen_limitation_ALAN_plant_decline_interaction_fitted":False,
    }
    dest=Path(out)
    require(not dest.exists(),"Refuse to overwrite empirical NPMS/GloPL join receipt")
    dest.mkdir(parents=True)
    taxon.to_csv(dest/"NPMS_recorded_species_x_GloPL_experiment_coverage.csv",
                 index=False)
    # Plot-year history is key to first GSHHG island/mainland and subsequent
    # PCNL comparison: export aggregated plot-level counts, not raw locations.
    plot=geo_sample.groupby("site_id",observed=True).agg(
        n_original_sample_visits=("id","nunique"),
        n_observed_years=("year","nunique"),
        first_sample_year=("year","min"),
        last_sample_year=("year","max"),
        lon=("lon","median"),lat=("lat","median"),
        one_source_country=("country","first")).reset_index()
    plot.to_csv(dest/"NPMS_true_survey_plot_repeats_and_location.csv",index=False)
    (dest/"NPMS_GloPL_source_exact_species_overlap_result.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_ACTUAL_NPMS_GLOPL_EXPERIMENTAL_POLLEN_LIMITATION_SPECIES_OVERLAP",
          json.dumps(report,sort_keys=True),flush=True)
    return report
def main():
    p=argparse.ArgumentParser()
    for name in ("package","glopl","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.package,a.glopl,a.out)
if __name__=="__main__":main()
