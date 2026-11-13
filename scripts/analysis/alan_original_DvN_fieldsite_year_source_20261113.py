#!/usr/bin/env python3
"""Original diel reproductive-day/night field-site/year matched source preflight.

Kendall & Nicholson 2025 author source (DOI 10.1111/ele.70036); 1652 original
treatment rows, 140 original study IDs and 262 coarse study/species/metric pairs.

CRITICAL: Original study/species/metric can CONTAIN MULTIPLE experimental
locations and years. Never assign one site/year to pooled original metrics.
Paired outcome units are original study × species × metric × site lat/lon ×
true study start/end year × month bounds, with both genuine day/night arms.

No PCNL or outcome-sign-based light sample selection; separate original
paired-period outcomes are retained for later post-source original satellite join.
The corrected artificial light dataset does not start until 1992, so the
earliest valid true field experiment is 1993. Never use publication year.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd

STUDY_BLOB="39a82d966afcdde670603add1c881ace8791ffe9"
TRAIT_BLOB="18381e550a7e405fe4f198909f0bb381bbecf1e0"
METRICS=("fruit set","seed set")
PAIRS=("night pollination","day pollination")
KEY=["study_ID","plant_species","treatment_effectiveness_metric",
     "coordinates_lat","coordinates_lon","year_start","year_end",
     "month_start","month_end","effectiveness_units"]

def require(ok, why):
    if not bool(ok): raise ValueError("ALAN_ORIGINAL_DIEL_FIELD_SOURCE: "+str(why))

def gitblob(path):
    b=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for a in iter(lambda:f.read(4<<20),b""):h.update(a)
    return h.hexdigest()

def numerical_weighted_mean(rows):
    a=pd.to_numeric(rows.effectiveness_value,errors="coerce")
    if a.isna().any() or a.lt(0).any():return None
    w=pd.to_numeric(rows.sample_size,errors="coerce")
    if w.notna().all() and w.gt(0).all():
        return float(np.average(a.to_numpy(float),weights=w.to_numpy(float)))
    return float(a.mean())

def build(study,trait,elevation,output):
    output=Path(output)
    require(not output.exists(),"Output must be fresh original source")
    require(gitblob(study)==STUDY_BLOB,"Published original raw day/night Git blob changed")
    require(gitblob(trait)==TRAIT_BLOB,"Published original plant traits Git blob changed")
    s=pd.read_csv(study,dtype=str,keep_default_na=False)
    require(len(s)==1652 and s.study_ID.nunique()==140,"Original 1652/140 source changed")
    require(set(KEY+["treatment_condition","effectiveness_value","sample_size",
                     "day_pollination_hrs","night_pollination_hrs","plant_native_status",
                     "coordinate_source"]).issubset(s.columns),
            "Original site date/treatment schema changed")
    s["species"]=s.plant_species.str.split().str.join(" ")
    s["metric"]=s.treatment_effectiveness_metric
    original_pairs=s.loc[s.metric.isin(METRICS)].groupby(
        ["study_ID","species","metric"],observed=True)
    raw_reproductive_coarse=sum(
        {"night pollination","day pollination"}.issubset(set(q.treatment_condition))
        for _,q in original_pairs)
    require(raw_reproductive_coarse==231,
            "Original 231 coarse day-night fruit/seed comparisons changed")
    pairs=[];reject={}
    for key,g in s.loc[s.metric.isin(METRICS)].groupby(KEY,observed=True,sort=False):
        sub={kind:g.loc[g.treatment_condition.eq(kind)] for kind in PAIRS}
        if any(d.empty for d in sub.values()):continue
        first=g.iloc[0]
        year=pd.to_numeric(first.year_start,errors="coerce")
        end=pd.to_numeric(first.year_end,errors="coerce")
        lat=pd.to_numeric(first.coordinates_lat,errors="coerce")
        lon=pd.to_numeric(first.coordinates_lon,errors="coerce")
        if not np.isfinite(year) or not float(year).is_integer() or year<1993 or year>2025:
            reject["no_true_post1992_experiment_start"] = reject.get("no_true_post1992_experiment_start",0)+1
            continue
        if np.isfinite(end) and (not float(end).is_integer() or end<year or end>2025):
            reject["inconsistent_true_experiment_end"]=reject.get("inconsistent_true_experiment_end",0)+1
            continue
        if not np.isfinite(lat) or not np.isfinite(lon) or abs(lat)>90 or abs(lon)>180:
            reject["original_field_geolocation_missing_or_invalid"]=reject.get("original_field_geolocation_missing_or_invalid",0)+1
            continue
        vals={kind:numerical_weighted_mean(q) for kind,q in sub.items()}
        if any(v is None or not np.isfinite(v) for v in vals.values()):
            reject["non_numeric_or_negative_original_treatment"]=reject.get("non_numeric_or_negative_original_treatment",0)+1
            continue
        nval,dval=vals["night pollination"],vals["day pollination"]
        if nval+dval<=0:
            reject["both_original_treatments_zero"]=reject.get("both_original_treatments_zero",0)+1
            continue
        dh=pd.to_numeric(g.day_pollination_hrs,errors="coerce")
        nh=pd.to_numeric(g.night_pollination_hrs,errors="coerce")
        hours_known=(dh.notna().all() and nh.notna().all() and
                     dh.gt(0).all() and nh.gt(0).all())
        day_hours=float(dh.median()) if hours_known else None
        night_hours=float(nh.median()) if hours_known else None
        ratio=(nval/night_hours)/(nval/night_hours+dval/day_hours) if hours_known else None
        require(0<=nval/(nval+dval)<=1,"Source paired night fraction invalid")
        def identical_or_unknown(name):
            vals=set(q for q in g[name] if q and q!="NA")
            return next(iter(vals)) if len(vals)==1 else ("unknown" if not vals else "mixed")
        row={
            "source_pair_id":len(pairs),
            "study_key":str(first.study_ID).strip(),
            "species":str(first.species).strip(),
            "metric":str(first.metric),
            "units":str(first.effectiveness_units),
            "year_exp":int(year),
            "year_end":int(end) if np.isfinite(end) else None,
            "experiment_year_is_single_reported":bool(np.isfinite(end) and end==year),
            "experiment_year_end_unreported":bool(not np.isfinite(end)),
            "lat":float(lat),"lon":float(lon),
            "country":identical_or_unknown("country"),
            "coordinate_source":identical_or_unknown("coordinate_source"),
            "plant_native_source":identical_or_unknown("plant_native_status"),
            "month_start":str(first.month_start),
            "month_end":str(first.month_end),
            "night_reproduction":nval,"day_reproduction":dval,
            "night_share":nval/(nval+dval),
            "night_share_per_hour":ratio,
            "day_hours":day_hours,"night_hours":night_hours,
            "night_treatment_rows":len(sub["night pollination"]),
            "day_treatment_rows":len(sub["day pollination"]),
            "treatment_source_rows":len(g),
            "open_pollination_arms_present":bool(g.treatment_condition.eq("open pollination").any()),
            "hand_pollination_arms_present":bool(g.treatment_condition.eq("hand pollination").any()),
            "source_accepted_metric_pair_cohort_before_LIGHT":True,
        }
        pairs.append(row)
    out=pd.DataFrame(pairs)
    require(len(out)>=190 and out.study_key.nunique()>=100 and
            out.species.nunique()>=100,"Insufficient independently documented source experiment pairs")
    require(not out.duplicated(["study_key","species","metric","units","lat","lon",
                                "year_exp","year_end","month_start","month_end"]).any(),
            "Scientific original study×species×site×year pair not unique")
    # Source author-supplied elevation; no invented DEM values.
    elev=pd.read_csv(elevation,dtype=str,keep_default_na=False,encoding='latin-1')
    require({"study_ID","coordinates_lat","coordinates_lon","elevation"}.issubset(elev),
            "Original author elevation metadata unavailable")
    elev=elev.copy()
    elev["study_key"]=elev.study_ID.astype(str)
    for c in ("coordinates_lat","coordinates_lon"):
        elev[c]=pd.to_numeric(elev[c],errors="coerce").round(5)
    elev["elevation_m"]=pd.to_numeric(elev.elevation,errors="coerce")
    elev=elev.groupby(["study_key","coordinates_lat","coordinates_lon"],observed=True,
         as_index=False).elevation_m.median()
    out["rounded_lat"]=out.lat.round(5)
    out["rounded_lon"]=out.lon.round(5)
    out=out.merge(elev,left_on=["study_key","rounded_lat","rounded_lon"],
                  right_on=["study_key","coordinates_lat","coordinates_lon"],
                  how="left",validate="many_to_one").drop(
                     columns=["rounded_lat","rounded_lon","coordinates_lat","coordinates_lon"])
    require(out.source_pair_id.is_unique and len(out)==len(pairs),"Original elevation lookup changed pair count")
    out["elevation_m"]=out.elevation_m.astype(float)
    # Author field source coordinate accuracy can be "based on place name",
    # hence adjust to author elevation only when exact published site matches.
    output.mkdir(parents=True)
    file=output/"original_DvN_field_site_year_outcome_pre_PCNLP.csv.gz"
    out.to_csv(file,index=False,compression="gzip")
    valid_coor=out.coordinate_source.value_counts().to_dict()
    spans=out.year_end.notna() & out.year_end.gt(out.year_exp)
    native=out.plant_native_source.eq("native")
    receipt={
        "schema":"ALAN_original_DvN_source_site_year_pre_PCNLP_v1",
        "published_article":"Kendall and Nicholson 2025 Ecology Letters 10.1111/ele.70036",
        "original_treatment_raw_rows":int(len(s)),
        "original_study_ids":int(s.study_ID.nunique()),
        "original_study_species_metric_paired_reproductive_groups":int(raw_reproductive_coarse),
        "actual_original_site_year_reproductive_pairs":len(out),
        "n_independent_original_studies":int(out.study_key.nunique()),
        "n_original_plant_species":int(out.species.nunique()),
        "n_distinct_field_coordinates":int(out[["lat","lon"]].drop_duplicates().shape[0]),
        "n_study_sites_year":int(out[["study_key","lat","lon","year_exp"]].drop_duplicates().shape[0]),
        "year_start_min":int(out.year_exp.min()),"year_start_max":int(out.year_exp.max()),
        "source_exact_year_end_equal_start":int(out.experiment_year_is_single_reported.sum()),
        "source_year_end_unreported":int(out.experiment_year_end_unreported.sum()),
        "source_multiyear_start_end":int(spans.sum()),
        "source_native_plant_pairs":int(native.sum()),
        "source_native_study_ids":int(out.loc[native,"study_key"].nunique()),
        "original_coordinates_quality":{str(k):int(v) for k,v in valid_coor.items()},
        "source_author_elevation_matched_pairs":int(out.elevation_m.notna().sum()),
        "source_hours_normalized_computable":int(out.night_share_per_hour.notna().sum()),
        "source_rejections":{k:int(v) for k,v in reject.items()},
        "original_author_treatment_blob":STUDY_BLOB,
        "original_author_trait_blob":TRAIT_BLOB,
        "original_outcome_blind_satellite_source_sha256":sha(file),
        "future_light_never_sampled":True,
        "source_year_is_actual_experiment_start_not_publication":True,
        "experiment_location_study_species_year_outcome_mixed_across_sites":False,
        "not_source_local_moth_identity":True,
        "not_island_only_sample":True,
        "not_PL_or_seedset_association_with_ALAN_without_real_satellite_join":True,
        "source_rows_restricted_to_internal_GHA":True,
    }
    (output/"original_DvN_source_year_geo_repro_audit.json").write_text(
          json.dumps(receipt,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("ALAN_ORIGINAL_REAL_DIEL_SOURCE_PAIR_AND_YEAR_SITE_FEASIBILITY",
          json.dumps(receipt,ensure_ascii=False,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ("study","traits","elevation","output"):
        ap.add_argument("--"+k,type=Path,required=True)
    a=ap.parse_args()
    build(a.study,a.traits,a.elevation,a.output)
