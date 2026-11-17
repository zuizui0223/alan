#!/usr/bin/env python3
"""Original EIDC NPMS 2015-22 habitat-change flag vs source-pinned ALAN panel.

The author field mutipleHabs (sic) is plot-wise: Y when more than one
broad habitat was reported between visits, N otherwise. This is an
observation-process and habitat-classification flag, NOT an experimentally
measured alteration in wetland water table nor a randomized light treatment.

Only exact matching original sample_id, location_id, survey_id and original
sampling date are accepted across the 2023 and 2026 published extracts.
No plant records or volunteer IDs are written in output, only aggregate counts.
"""
import argparse
import json
import zipfile
from pathlib import Path

import pandas as pd

from scripts.acquisition import alan_npms_original_habitat_source_gate_20261117 as npms
from scripts.analysis import alan_npms_true_wetland_moisture_priorlight_20261117 as wet
from scripts.analysis import alan_npms_wetland_light_spatial_falsification_20261117 as check

SCHEMA="alan_npms_official_2015_22_habitat_change_panel_gate_v1"
DOI="10.5285/a12b5c70-c0e2-48cf-b0b7-3e95791feca9"
MEMBER="data/npmshabitatsamples_2015to2022.csv"
CRC32="16b150ac"
SIZE=4389198

def require(cond,msg):
    if not bool(cond):raise RuntimeError("ALAN_NPMS_HAB_CHANGE_STOP: "+msg)

def source_compare(habitat_zip,npms_zip):
    with zipfile.ZipFile(habitat_zip) as z:
        require(z.testzip() is None,"original habitat release failed CRC test")
        info=z.getinfo(MEMBER)
        require(info.file_size==SIZE and f"{info.CRC:08x}"==CRC32,
                "original published NPMS habitat member changed")
        with z.open(MEMBER) as f:
            old=pd.read_csv(f,dtype=str,keep_default_na=False)
    require(len(old)==21185 and old.location_id.nunique()==4712,
            "actual 2015-22 source cardinality changed")
    require(set(old.mutipleHabs)=={"N","Y"},"source flag changed from Y/N")
    old=old.drop_duplicates()
    require(len(old)==21184 and old.sample_id.is_unique,"original duplicate gate changed")
    old["date"]=pd.to_datetime(old.date_start,utc=True,errors="coerce")
    require(old.date.notna().all() and old.date.dt.year.between(2015,2022).all(),
            "invalid source survey years")
    by=old.groupby("location_id",sort=False).agg(
        codes=("mutipleHabs","nunique"),flag=("mutipleHabs","first"),
        distinct_broad=("NPMS_broad_habitat",lambda x:x[x.ne("")].nunique()),
        visited_years=("date",lambda x:x.dt.year.nunique()))
    require(len(by)==4712 and (by.codes==1).all(),"source site flag inconsistency")
    require(int((by.flag=="Y").sum())==351,"published changed plots drift")
    require(((by.flag=="Y")==(by.distinct_broad>=2)).all(),
            "author change flag does not coincide with two reported habitats")
    with zipfile.ZipFile(npms_zip) as z:
        require(z.testzip() is None,"NPMS 2025 original archive CRC")
        item=z.getinfo("data/sampleinfowithlatlong_2015to2025.csv")
        require(item.file_size==4095902 and f"{item.CRC:08x}"=="5e52f384",
                "2025 source exact sample CSV changed")
        with z.open(item) as f:now=pd.read_csv(f,dtype=str,keep_default_na=False)
    require(len(now)==27611 and now.id.is_unique,
            "NPMS original 2025 survey keys changed")
    joint=old.merge(now[["id","location_id","survey_id","date_start"]],
         left_on="sample_id",right_on="id",how="left",indicator=True,
         suffixes=("_old","_new"),validate="one_to_one")
    old_date=pd.to_datetime(joint.date_start_old,utc=True,errors="coerce")
    new_date=pd.to_datetime(joint.date_start_new,utc=True,errors="coerce")
    passrow=(joint["_merge"].eq("both")&
        joint.location_id_old.eq(joint.location_id_new)&
        joint.survey_id_old.eq(joint.survey_id_new)&old_date.eq(new_date))
    good_sites=set(joint.loc[passrow,"location_id_old"])
    bad_sites=set(joint.loc[~passrow,"location_id_old"])
    coherent=good_sites-bad_sites
    by["version_ok"]=by.index.isin(coherent)
    meta={
      "schema":SCHEMA,"EIDC_habitat_change_DOI":DOI,
      "habitat_member_crc32":CRC32,"habitat_member_original_bytes":SIZE,
      "actual_source_rows":21185,"distinct_source_samples":21184,
      "physical_plots_all":4712,
      "physical_plots_ever_source_reported_broad_change":351,
      "physical_plots_only_one_source_reported_broad_category":4361,
      "cross_version_sample_rows_identical":int(passrow.sum()),
      "cross_version_sample_rows_unmatched_or_inconsistent":int((~passrow).sum()),
      "cross_version_physical_plots_consistent":len(coherent),
      "cross_version_physical_plots_rejected_for_any_record_mismatch":len(bad_sites),
      "source_version_is_observer_classification_not_new_field_experiment":True,
      "missing_or_nonreported_species_not_assumed_absent":True,
    }
    return by,meta

def selected_support(common,by):
    visits=common[["location_id","inland_wetland"]].drop_duplicates()
    require(visits.location_id.is_unique,"source wetland identity inconsistent")
    info=visits.merge(by,left_on="location_id",right_index=True,
                     how="left",validate="one_to_one")
    out={"total_source_matched_physical_plots":len(visits)}
    for marker,g in info.groupby("inland_wetland"):
        z=g.loc[g.version_ok.eq(True)]
        out["wetland" if marker else "other"]={
            "original_fixed_plots":len(g),
            "exact_crossversion_ID_plots":len(z),
            "not_exactly_linked_or_missing_plots":len(g)-len(z),
            "source_reported_change_Y_plots":int(z.flag.eq("Y").sum()),
            "source_reported_no_change_N_plots":int(z.flag.eq("N").sum()),
            "no_change_with_2plus_source_years":int(
                (z.flag.eq("N") & z.visited_years.ge(2)).sum())
        }
    return out

def model(common,valid_sites,name):
    sample=common.loc[common.location_id.isin(valid_sites)].copy()
    sample=check.mixed_reselect(sample)
    res=check.estimate(sample,name)
    point=res.get("triple_habitat_moisture_light",{})
    return {
        "status":res["status"],
        "plot_support":res["support"],
        "beta":point.get("beta"),
        "plot_cluster_95CI":point.get("plot_95CI"),
        "two_way":point.get("two_way")
    }

def execute(habitat_2022,npms_2025,cover,receipt,eiv):
    by,metadata=source_compare(habitat_2022,npms_2025)
    src,common,eligible=wet.study(
        npms_2025,cover,receipt,eiv,return_source_panels=True)
    require(len(common)==11902 and common.location_id.nunique()==326,
            "original positive-cover wetland study frame changed")
    out={**metadata,
         "original_NPMS_positive_cover_rows":len(common),
         "source_2015_22_habitat_change_support":selected_support(common,by),
         "posthoc_source_version_sensitivities":{},
         "interpretation":"source observation-habitat change not wetland drainage or ALAN causation"
    }
    cases={
        "only_exact_2023_vs_2026_matched_plot_IDs":
            by.index[by.version_ok],
        "exclude_author_reported_broad_habitat_change":
            by.index[by.version_ok & by.flag.eq("N")],
        "exclude_change_require_two_source_survey_years":
            by.index[by.version_ok & by.flag.eq("N") & by.visited_years.ge(2)],
        "flagged_broad_habitat_change_source_support_only":
            by.index[by.version_ok & by.flag.eq("Y")]
    }
    for name,subset in cases.items():
        out["posthoc_source_version_sensitivities"][name]=model(
            common,set(subset),name)

    # Joint source restrictions: no author-reported broad habitat change,
    # *and* unchanged original botanical status across the same plot/plant/
    # survey-method/season panel. Both are retrospective audit restrictions.
    baseline=set(cases["exclude_change_require_two_source_survey_years"])
    official_N=common.loc[common.location_id.isin(baseline)].copy()
    source_status=official_N.groupby("panel_key").source_record_status.nunique(dropna=False)
    stable_keys=set(source_status.index[source_status==1])
    jointly=official_N.loc[official_N.panel_key.isin(stable_keys)].copy()
    jointly=check.mixed_reselect(jointly)
    out["posthoc_source_version_sensitivities"][
        "official_unchanged_and_original_record_status_constant"
    ]=model(jointly,set(jointly.location_id),
            "official_unchanged_and_original_record_status_constant")
    out["joint_source_record_status_gates"]={
        "original_plot_species_survey_season_panels_with_known_official_stable_habitat":
            len(source_status),
        "panels_where_source_original_record_status_constant":len(stable_keys),
        "panels_where_source_original_record_status_changed":
            int((source_status>1).sum())
    }

    # The independent official change release ends at 2022. Restrict the
    # botanical outcome window as a sensitivity to that exact time horizon.
    strict_2015_22=official_N.loc[official_N.year.le(2022)].copy()
    strict_2015_22=check.mixed_reselect(strict_2015_22)
    out["posthoc_source_version_sensitivities"][
        "official_unchanged_and_outcome_year_le2022"
    ]=model(strict_2015_22,set(strict_2015_22.location_id),
            "official_unchanged_and_outcome_year_le2022")

    # Species-composition countercheck within the official stable-source
    # matched sites, not between differently drawn case/control species.
    aa=set(official_N.loc[official_N.inland_wetland,"species"])
    bb=set(official_N.loc[~official_N.inland_wetland,"species"])
    same_species=aa&bb
    shared=official_N.loc[official_N.species.isin(same_species)].copy()
    shared=check.mixed_reselect(shared)
    out["source_stable_shared_species_count"]=len(same_species)
    out["posthoc_source_version_sensitivities"][
        "official_unchanged_same_taxa_in_both_habitats"
    ]=model(shared,set(shared.location_id),
            "official_unchanged_same_taxa_in_both_habitats")

    # Stronger habitat-specific temporal control. A species may change
    # differently in a wetland and a dry habitat in a given calendar year.
    # This replaces the old species-year FE with species-by-habitat-year.
    # The resulting coefficient is a falsification test of the interpretation.
    official_N=check.mixed_reselect(official_N)
    official_N["species_habitat_year"]=(
        official_N.species+"|"+official_N.inland_wetland.astype(int).astype(str)
        +"|"+official_N.year.astype(str))
    base=check.common_support(official_N)
    temporal={"source_support":base}
    if base["pass"]:
        from scripts.analysis import alan_npms_independent_moisture_priorlight_20261117 as core
        fitted=core.fit(official_N,["light_moist","light_moist_wetland","doy30"],
             ["panel_key","species_habitat_year","plotyear_method_season"],
             two_way=True)
        temporal["model_status"]=fitted["status"]
        if fitted["status"]=="observational_positive_cover_only":
            term=fitted["coefficients"]["light_moist_wetland"]
            temporal.update(beta=term["beta"],nominal_plot_95CI=term["nominal_95CI"],
              plot_p=term["nominal_two_sided_p"],
              two_way=fitted["two_way_plot_species_robustness"])
        else:
            temporal["stop_reason"]="original outcome has no remaining within-source variation"
    else:
        temporal["model_status"]="original_wetland_common_support_STOP"
    out["posthoc_falsifier_species_by_habitat_by_year_FE"]=temporal
    out["no_causal_hydrology_or_light_mediation_estimated"]=True
    out["reporting_restrictions"]=[
        "Version-joined sample data from separate 2023/2026 NPMS releases, not distinct field sampling programmes",
        "Habitat change Y is a retrospective site-level category change across source visits 2015-2022",
        "N means no *reported* broad category change within available visits, not actual habitat stability",
        "Minimum independent wetland site common-support gates block unsupported new regression",
        "Original plant cover was only measured conditional on a positive species record",
        "The stage does not test pollen limitation, plant disappearance, flower visitation or manipulated lighting"
    ]
    return out

def main():
    a=argparse.ArgumentParser()
    for x in ("habitat-2022","npms-2025","cover","receipt","eiv","out"):
        a.add_argument("--"+x,required=True)
    q=a.parse_args()
    r=execute(q.habitat_2022,q.npms_2025,q.cover,q.receipt,q.eiv)
    f=Path(q.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(r,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("NPMS_HABITAT_CHANGE_EXACT_ID_FALSIFICATION",json.dumps({
      "source":{k:r[k] for k in ("physical_plots_all",
        "physical_plots_ever_source_reported_broad_change",
        "cross_version_sample_rows_identical",
        "cross_version_physical_plots_consistent")},
      "cohort":r["source_2015_22_habitat_change_support"],
      "models":{k:{"status":v["status"],
                    "wet_sites":v["plot_support"]["inland_wetland"]["plots"],
                    "beta":v["beta"],"CI":v["plot_cluster_95CI"]}
                for k,v in r["posthoc_source_version_sensitivities"].items()}
    },sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
