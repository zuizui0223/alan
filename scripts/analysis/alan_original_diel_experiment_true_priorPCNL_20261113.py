#!/usr/bin/env python3
"""Actual within-study day/night plant reproduction experiment × previous annual PCNL.

Primary: same original study_ID + study plant binomial + benefit metric, night
pollination versus day pollination exclusion arms. Retain only exact shared
experiment year, original coordinates and measurement unit, 1993–2022.
Independent previous-calendar-year Li–Cao corrected annual PCNL 1992–2021,
MD5 verified for every actual annual raster, sampled at source study site.

This is a cross-study observational comparison of directly measured source
day/night reproductive responses and prior satellite night radiance, NOT an
experimental intervention on light, proof of moth-only service, or plant
community change. The authored exposure and study site metadata remain
ephemeral; publish only aggregate results.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import t

from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    fetch_one, manifest_read
)

SCHEMA="ALAN_true_diel_daynight_exclusion_reproductive_field_priorPCNL_1993_2022_v1"
AUTHOR_BLOB="39a82d966afcdde670603add1c881ace8791ffe9"
TRAITS_BLOB="18381e550a7e405fe4f198909f0bb381bbecf1e0"
YEARS=tuple(range(1992,2022))
SHARDS=6
PAIR_N=146
STUDY_N=87
SPECIES_N=94
REPRO=("fruit set","seed set")
EXPOSURE="pcnl_prior_mean_log1p"
KEY=["study","species","metric"]


def require(cond,why):
    if not cond:raise ValueError("ALAN_DIEL_REPRO_PCNl: "+str(why))


def gitblob(path):
    b=Path(path).read_bytes()
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as fp:
        for block in iter(lambda: fp.read(4<<20),b""):
            h.update(block)
    return h.hexdigest()


def author_pairs(study_path,traits_path):
    """Source-original exact within-study paired outcomes, not species-level
    external day/night labels or GloBI pollination guesses."""
    require(gitblob(study_path)==AUTHOR_BLOB,"Original study treatment Git blob drift")
    require(gitblob(traits_path)==TRAITS_BLOB,"Original study botanical trait blob drift")
    raw=pd.read_csv(study_path,dtype=str,keep_default_na=False)
    tr=pd.read_csv(traits_path,dtype=str,keep_default_na=False)
    require(len(raw)==1652 and len(tr)==142,"Original publication rows drift")
    cols={"study_ID","plant_species","treatment_effectiveness_metric",
          "treatment_condition","effectiveness_value","effectiveness_units",
          "year_start","year_end","coordinates_lat","coordinates_lon",
          "sample_size","plant_native_status","coordinate_source",
          "day_pollination_hrs","night_pollination_hrs"}
    require(cols.issubset(raw.columns),"Required original study metadata absent")
    original=0
    docs=[]
    for (study,species,metric),g in raw.groupby(
            ["study_ID","plant_species","treatment_effectiveness_metric"],
            observed=True,sort=True):
        day=g.loc[g.treatment_condition.eq("day pollination")]
        night=g.loc[g.treatment_condition.eq("night pollination")]
        if day.empty or night.empty or metric not in REPRO:continue
        original+=1
        pair=pd.concat([day,night],ignore_index=True)
        # Confirm the EXACT same published study/population/year and
        # original response scale; a paper's several sites/years cannot be
        # silently averaged as an exposure at an invented reference site.
        if (pair.year_start.nunique()!=1 or pair.year_end.nunique()!=1 or
            pair.coordinates_lat.nunique()!=1 or
            pair.coordinates_lon.nunique()!=1 or
            pair.effectiveness_units.nunique()!=1 or
            pair.exclude.nunique()!=1 or pair.exclude.iloc[0]!="n"):
            continue
        sy=pair.year_start.iloc[0]
        ey=pair.year_end.iloc[0]
        if not sy.isdigit() or len(sy)!=4 or not (1993<=int(sy)<=2022):
            continue
        if ey not in ("NA","",sy):continue
        latitude=pd.to_numeric(pair.coordinates_lat.iloc[0],errors="coerce")
        longitude=pd.to_numeric(pair.coordinates_lon.iloc[0],errors="coerce")
        if not (np.isfinite(latitude) and np.isfinite(longitude) and
                -90<=latitude<=90 and -180<=longitude<=180):
            continue
        def observed_mean(d):
            v=pd.to_numeric(d.effectiveness_value,errors="coerce")
            if not np.isfinite(v).all() or (v<0).any():return np.nan
            sample=pd.to_numeric(d.sample_size,errors="coerce")
            if (np.isfinite(sample).all() and (sample>0).all()):
                return float(np.average(v,weights=sample))
            return float(v.mean())
        dval=observed_mean(day);nval=observed_mean(night)
        if not (np.isfinite(dval) and np.isfinite(nval) and
                dval>=0 and nval>=0 and (dval+nval)>0):
            continue
        # Exact day+night clock hours are a secondary model; do not impute
        # missing hours or assume seed output is linear in pollination hours.
        dhr=pd.to_numeric(pair.day_pollination_hrs,errors="coerce")
        nhr=pd.to_numeric(pair.night_pollination_hrs,errors="coerce")
        if (np.isfinite(dhr).all() and np.isfinite(nhr).all() and
            (dhr>0).all() and (nhr>0).all() and
            np.ptp(dhr.to_numpy())<1e-7 and np.ptp(nhr.to_numpy())<1e-7):
            hours_day=float(dhr.iloc[0]);hours_night=float(nhr.iloc[0])
        else:
            hours_day=np.nan;hours_night=np.nan
        nightshare=float(nval/(nval+dval))
        original_native=pair.plant_native_status.unique()
        source_coord=pair.coordinate_source.unique()
        docs.append({
           "study":study,"species":species,"metric":metric,
           "year_exp":int(sy),"lat":float(latitude),"lon":float(longitude),
           "original_unit":pair.effectiveness_units.iloc[0],
           "native_status":original_native[0] if len(original_native)==1 else "unknown",
           "coord_quality":source_coord[0] if len(source_coord)==1 else "mixed",
           "day_value":dval,"night_value":nval,
           "day_hours":hours_day,"night_hours":hours_night,
           "night_share":nightshare,
           "night_minus_day":2*nightshare-1,
           "night_per_hour_share": (
                 (nval/hours_night)/(nval/hours_night+dval/hours_day)
                 if np.isfinite(hours_day) and np.isfinite(hours_night) else np.nan),
           "original_day_rows":len(day),"original_night_rows":len(night)
        })
    require(original==231,"Original author source reproductive day/night pair count drift")
    q=pd.DataFrame(docs).sort_values(KEY).reset_index(drop=True)
    require(len(q)==PAIR_N and q.study.nunique()==STUDY_N and
            q.species.nunique()==SPECIES_N,
            "Original strict same-year/same-site author study population drift")
    require(q[KEY].is_unique if hasattr(q[KEY],"is_unique") else not q.duplicated(KEY).any(),
            "Duplicate published day/night comparison group")
    q["row_id"]=np.arange(len(q))
    q["has_hours"]=q.night_per_hour_share.notna()
    q["lat2"]=q.lat**2
    q["lon2"]=q.lon**2
    q["geo10"]=((q.lat+90)//10).astype(int).astype(str)+"_"+(
               ((q.lon+180)//10).astype(int).astype(str))
    require(q.year_exp.between(1993,2022).all(),"Source experiment year invalid")
    require(q[["lat","lon","night_share"]].notna().all().all(),"Source exact metric missing")

    tr["species"]=tr.plant_species.astype(str).str.strip()
    tr=tr[["species","bloom_period_simple","breeding_system"]].copy()
    require(tr.species.is_unique,"Author traits duplicate original botanical name")
    q=q.merge(tr,on="species",how="left",validate="many_to_one")
    return q,{
       "original_author_study_treatment_rows":len(raw),
       "original_author_reproductive_paired_groups":original,
       "strict_author_same_study_species_metric_year_site_unit_groups":len(q),
       "strict_author_distinct_studies":int(q.study.nunique()),
       "strict_author_distinct_species":int(q.species.nunique()),
       "strict_native_source_groups":int(q.native_status.eq("native").sum()),
       "strict_nonnative_source_groups":int(q.native_status.isin(
            ("non-native","nonnative")).sum()),
       "strict_source_exact_years":[int(q.year_exp.min()),int(q.year_exp.max())],
       "reported_original_experiment_hours_pairs":int(q.has_hours.sum()),
       "study_pairs_per_publication_median":float(q.groupby("study").size().median()),
       "studies_with_2plus_source_comparisons":int(q.groupby("study").size().ge(2).sum()),
       "source_original_author_published_study_gitblob":AUTHOR_BLOB,
       "source_original_author_published_trait_gitblob":TRAITS_BLOB,
       "original_study_year_not_PUBLICATION_YEAR":True,
       "same_study_species_reproduction_metric_scale_matched":True,
       "no_abstract_moth_identity_assigned":True
    }


def sample_one(tif,x):
    """Original published PCNL pixel (EPSG native reprojection); no
    inverse-distance fill, zero luminance still valid, masked stays NA."""
    from rasterio.warp import transform
    xs=x.lon.to_numpy(float);ys=x.lat.to_numpy(float)
    if tif.crs.to_epsg()!=4326:
        xs,ys=transform("EPSG:4326",tif.crs,xs.tolist(),ys.tolist())
    coords=list(zip(xs,ys))
    r=tif.sample(coords,indexes=1,masked=True)
    values=np.array([np.nan if np.ma.is_masked(t[0]) else float(t[0])
                     for t in r],float)
    return np.where(np.isfinite(values)&(values>=0),values,np.nan)


def shard_years(shard):
    require(shard in range(SHARDS),"Not an original source year shard")
    return [y for y in YEARS if (y-1992)%SHARDS==shard]


def worker(study,traits,manifest,shard,folder,out):
    import rasterio
    x,meta=author_pairs(study,traits)
    files=manifest_read(manifest)["files"]
    require(not Path(out).exists(),"Source partial output already present")
    Path(out).mkdir(parents=True)
    Path(folder).mkdir(parents=True,exist_ok=True)
    sums=np.zeros(len(x),float)
    nobs=np.zeros(len(x),np.int16)
    receipts=[]
    for year in shard_years(shard):
        file="PCNL"+str(year)+".tif"
        source=fetch_one(file,files[file],Path(folder)/file)
        with rasterio.open(Path(folder)/file) as tif:
            require(tif.count==1 and tif.crs is not None,
                    "Original annual PCNL not a valid georeferenced single band")
            px=sample_one(tif,x)
        eligible=x.year_exp.to_numpy(int)>year
        good=eligible&np.isfinite(px)
        sums[good]+=np.log1p(px[good]);nobs[good]+=1
        receipts.append({"year":year,
                         "n_source_pairs_before_year":int(eligible.sum()),
                         "n_valid_eligible_original_pixels":int(good.sum()),
                         **source})
        print("ORIGINAL_DIEL_PCNL_YEAR_CHECK",json.dumps(receipts[-1]),flush=True)
    dest=Path(out)/("diel_pcNL_shard_"+str(shard)+".csv.gz")
    pd.DataFrame({"row_id":x.row_id,
                  "partial_n":nobs,
                  "partial_sum_log1p_PCNL":sums}).to_csv(
                      dest,index=False,compression="gzip")
    receipt={"schema":SCHEMA,"shard":shard,
             "original_author_gitblob":AUTHOR_BLOB,
             "original_author_trait_gitblob":TRAITS_BLOB,
             "strict_pair_rows":PAIR_N,
             "verified_years":shard_years(shard),
             "original_md5_verified_years":receipts,
             "partial_SHA256":sha256(dest),
             "no_future_PCNL_for_author_experiment":True}
    (Path(out)/("diel_pcNL_shard_"+str(shard)+".json")).write_text(
         json.dumps(receipt,indent=2,sort_keys=True)+"\n")


def wls_study(g,formula,*,term=EXPOSURE):
    """Clusters by original independent published study; no pseudo
    replicate CI for multiple metrics from one botanical study."""
    meta={"n_pairs":len(g),"n_studies":int(g.study.nunique()),
          "n_species":int(g.species.nunique()),
          "n_source_light_regions_10deg":int(g.geo10.nunique()),
          "n_publications_2plus_pairs":int(g.groupby("study").size().ge(2).sum())}
    if len(g)<60 or g.study.nunique()<45 or g[EXPOSURE].std()<.03:
        return {"status":"NOT_ESTIMABLE_STUDY_LIGHT_SUPPORT",**meta}
    q=g.copy()
    q["study_equal_weight"]=1/q.groupby("study").study.transform("size")
    m=smf.wls(formula,data=q,weights=q.study_equal_weight).fit(
          cov_type="cluster",cov_kwds={"groups":q.study},use_t=False)
    rank=np.linalg.matrix_rank(m.model.exog)
    if (rank<len(m.params) or term not in m.params or
        not np.isfinite(float(m.bse[term]))):
        return {"status":"NOT_ESTIMABLE_FE_DESIGN_RANK",
                **meta,"rank":int(rank),"n_parameters":len(m.params)}
    ci=m.conf_int().loc[term]
    return {"status":"EXPLORATORY_WITHIN_STUDY_PAIRED_REPRODUCTION_VS_PRIOR_PCNL",
            "beta_nightshare_pp_per_mean_logPCNL":float(m.params[term]),
            "95CI_publication_cluster":[float(ci.iloc[0]),float(ci.iloc[1])],
            "p_nominal_two_sided":float(m.pvalues[term]),
            "R2":float(m.rsquared),"formula":formula,
            "n_parameters":len(m.params),
            **meta,"response":"100*night/(night+day) same study species field and metric",
            "light_observational_not_manipulated":True}


def fitted_models(x):
    x=x.copy()
    x["nightshare_pp"]=100*x.night_share
    x["hours_pp"]=100*x.night_per_hour_share
    x["log_reported_night_vs_day_hours"]=np.log(
        x.night_hours/x.day_hours)
    x["source_year"]=x.year_exp
    x["coord_direct"]=x.coord_quality.eq("provided in text")
    for k in ["bloom_period_simple","native_status"]:
        x[k]=x[k].fillna("unknown").replace({"NA":"unknown","":"unknown"})
    x["source_20deg"]=((x.lat+90)//20).astype(int).astype(str)+"_"+(
                       ((x.lon+180)//20).astype(int).astype(str))
    a="nightshare_pp ~ "+EXPOSURE+" + source_year + lat + lat2 + lon + lon2 + C(metric)"
    result={}
    variants={
        "primary_original_paired_study_broad_geography":(x,a),
        "20deg_georegion_FE":(x,a+" + C(source_20deg)"),
        "original_recorded_native_only":(x.loc[x.native_status.eq("native")],a),
        "location_directly_stated_in_paper_only":(x.loc[x.coord_direct],a),
        "one_source_metric_fruit_set_only":(x.loc[x.metric.eq("fruit set")],a),
        "one_source_metric_seed_set_only":(x.loc[x.metric.eq("seed set")],a),
        "hours_available_no_hours_normalization":(
            x.loc[x.night_per_hour_share.notna()],a+
            " + log_reported_night_vs_day_hours"),
        "night_share_PER_REPORTED_HOUR_sensitivity":(
            x.loc[x.night_per_hour_share.notna()],
            a.replace("nightshare_pp","hours_pp")),
        "reported_night_blooming_author_trait_only":(
            x.loc[x.bloom_period_simple.eq("night")],a)
    }
    for name,(g,formula) in variants.items():
        print("ALAN_ACTUAL_DIEL_PCNL_FIT",name,len(g),flush=True)
        try:result[name]=wls_study(g,formula)
        except (ValueError,np.linalg.LinAlgError) as e:
            result[name]={"status":"NOT_ESTIMABLE_MODEL_FAILURE","reason":str(e)[:220]}
    # A genuine WITHIN original-study light regression is only possible
    # if papers contain multiple field sites or experiment years AND both
    # original year, exposure and actual treatment outcomes are observed.
    studylight=x.groupby("study").agg(
       n_pairs=("row_id","size"),
       n_coords=("lat","nunique"),
       source_light_min=(EXPOSURE,"min"),
       source_light_max=(EXPOSURE,"max"))
    variable=studylight.loc[(studylight.source_light_max-
                             studylight.source_light_min).gt(.03)]
    src=x[x.study.isin(variable.index)]
    result["actual_within_published_study_source_support"]={
       "studies_with_2plus_genuine_source_pairs":int(studylight.n_pairs.ge(2).sum()),
       "studies_with_original_historical_light_span_gt_003":int(len(variable)),
       "cases_in_these_studies":int(len(src)),
       "distinct_source_species_in_these_studies":int(src.species.nunique()),
       "distinct_geographic_study_locations_exposure_variation_not_identical_light":True
    }
    if len(src)>=75 and src.study.nunique()>=25:
        result["within_publication_fixed_effects"]=wls_study(
           src,a+" + C(study)")
    else:
        result["within_publication_fixed_effects"]={
           "status":"NOT_ESTIMABLE_INDEPENDENT_STUDY_LIGHT_SUPPORT",
           "study_support":result["actual_within_published_study_source_support"]}
    result["no_actual_moth_identity_or_ALAN_intervention_in_author_study"]=True
    return result


def collect(study,traits,parts,out):
    if Path(out).exists():raise FileExistsError(out)
    x,meta=author_pairs(study,traits)
    n=np.zeros(len(x),np.int16)
    sums=np.zeros(len(x),float)
    receipts=[]
    years=[]
    for i in range(SHARDS):
        f=Path(parts)/("diel_pcNL_shard_"+str(i)+".csv.gz")
        p=Path(parts)/("diel_pcNL_shard_"+str(i)+".json")
        require(f.is_file() and p.is_file(),"Original annual source shard missing")
        r=json.loads(p.read_text())
        require(r["schema"]==SCHEMA and r["shard"]==i and
                r["strict_pair_rows"]==PAIR_N and
                r["original_author_gitblob"]==AUTHOR_BLOB and
                r["original_author_trait_gitblob"]==TRAITS_BLOB and
                r["verified_years"]==shard_years(i) and
                sha256(f)==r["partial_SHA256"],
                "Original author study or satellite shard changed")
        d=pd.read_csv(f)
        require(d.row_id.eq(x.row_id).all() and len(d)==len(x),
                "Original source paired result ordering inconsistent")
        n+=d.partial_n.to_numpy(np.int16)
        sums+=d.partial_sum_log1p_PCNL.to_numpy(float)
        years.extend(r["verified_years"])
        receipts.append(r)
    require(sorted(years)==list(YEARS) and len(set(years))==len(YEARS),
            "Original 1992–2021 preceding-nightlight source incomplete")
    require(np.all(n<=x.year_exp.to_numpy()-1992),"Future year sampled for experimental study")
    complete=n==x.year_exp.to_numpy()-1992
    x[EXPOSURE]=np.where(complete,sums/np.maximum(1,x.year_exp-1992),np.nan)
    sample=x.loc[complete].copy()
    support={
       "original_source_pairs":PAIR_N,
       "original_studies":STUDY_N,
       "original_species":SPECIES_N,
       "source_years_verified":len(YEARS),
       "first_year":min(YEARS),"last_year":max(YEARS),
       "source_original_calibrated_raster_count":len(YEARS),
       "actual_valid_complete_PCNL_reproductive_comparisons":int(len(sample)),
       "actual_valid_original_studies":int(sample.study.nunique()),
       "actual_valid_original_species":int(sample.species.nunique()),
       "valid_native_study_comparisons":int(sample.native_status.eq("native").sum()),
       "valid_reported_daynight_hours_comparisons":int(sample.night_per_hour_share.notna().sum()),
       "valid_original_coordinate_direct_study_comparisons":int(sample.coord_quality.eq("provided in text").sum()),
       "actual_light_zero_cases":int(sample[EXPOSURE].eq(0).sum()),
       "actual_light_nonzero_cases":int(sample[EXPOSURE].gt(0).sum()),
       "actually_distinct_field_study_sites":int(sample[["lat","lon"]].drop_duplicates().shape[0]),
       "actual_geo10_source_regions":int(sample.geo10.nunique()),
       "original_study_src_and_repro":meta
    }
    models=fitted_models(sample) if len(sample)>0 else {
           "status":"NOT_ESTIMABLE_NO_GENUINE_PRIOR_LIGHT_PIXELS"}
    result={
        "schema":SCHEMA,"source_support":support,
        "models":models,
        "light_source":"Li–Cao Zenodo 10.5281/zenodo.17013414 original 1992-2021 PCNL annual MD5-verified calibrated satellite pixel mean log1p STRICTLY BEFORE author original true field experiment year",
        "author_study_original_Git_blob":AUTHOR_BLOB,
        "author_botanical_traits_original_Git_blob":TRAITS_BLOB,
        "no_causal_inference":{
            "within_source_experiment_day_night_reproduction_measured":True,
            "no_field_experiment_randomizes_artificial_nightlight":True,
            "source_pair_source_year_coordinate_matched":True,
            "no_moth_specific_pollinator_taxa_identified":True,
            "no_local_longitudinal_reproductive_change_observed":True,
            "author_published_native_status_not_wild_native_island_status":True,
            "many_studies_have_single_original_location_light":True,
            "sensitivity_analysis_not_independent_replication":True,
            "no_field_population_trait_by_location_exported":True
        }
    }
    Path(out).parent.mkdir(parents=True,exist_ok=True)
    Path(out).write_text(json.dumps(result,sort_keys=True,indent=2,
                                    ensure_ascii=False,allow_nan=False)+"\n")
    summary={k:{"status":v.get("status"),"beta":v.get(
        "beta_nightshare_pp_per_mean_logPCNL"),
        "CI":v.get("95CI_publication_cluster"),
        "p":v.get("p_nominal_two_sided"),
        "n":v.get("n_pairs"),"studies":v.get("n_studies")}
       for k,v in models.items() if isinstance(v,dict)}
    print("ALAN_ACTUAL_PAIRED_DIEL_REPRODUCTIVE_PCNL_RESULT",
          json.dumps({"source_support":support,"models":summary},
                     ensure_ascii=False),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="phase",required=True)
    a=sub.add_parser("worker")
    for k in ("study","traits","manifest","raster-dir","out"):
        a.add_argument("--"+k,type=Path,required=True)
    a.add_argument("--shard",type=int,required=True)
    b=sub.add_parser("collect")
    for k in ("study","traits","parts","out"):
        b.add_argument("--"+k,type=Path,required=True)
    z=p.parse_args()
    if z.phase=="worker":
        worker(z.study,z.traits,z.manifest,z.shard,
               getattr(z,"raster_dir"),z.out)
    else:
        collect(z.study,z.traits,z.parts,z.out)

if __name__=="__main__":main()
