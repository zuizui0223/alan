#!/usr/bin/env python3
"""Original 6-lighting-zone Xue 2026 plant species functional divergence.

Original ScienceDB 10.57760/sciencedb.35715, Monoculture_Growth.xlsx
original SHA256: 8a18f7a30fe8e4656495db511f4e19b0f7e97657233111ba175c7c4007ed4f39.
Original source-extracted 77 surviving cylinders, 6 physically distinct ALAN
lighting zones (3 lit, 3 dark), 4 original monoculture macrophyte taxa.
Original 24 *zone x species* means are pinned to exact Git blob SHA1.
No 77 independent light replicates, no 4 independent experiments.
Plants are cultivated in separate monoculture pots: relative productivity
is NOT observed competition. Survivor selection is not zero-filled.

Question: is the light-vs-dark change in source unit-level relative biomass
greater for invasive Myriophyllum aquaticum than for the native plants?
For each real lighting zone, I_z =
    log(invasive biomass) - average(log(each 3 native species biomass)).
The difference in mean I_z lit minus dark removes large between-species
baseline biomass scales, while preserving real within-zone cross-species
covariance; it is not a mixed-community invasive fraction.

Exact 3-of-6 *zone label exchangeability* null (20 total allocations).
No claim of randomized light-zone assignment without original allocation
ledger. Two-sided exact minimum is .10, regardless of impressive effects.
4999 stratified (3+3 physical lighting zone) bootstrap interval is
descriptive conditional on zones, not across experiments/countries.
All no-snail/single-species source biomass effects are distinct from the
multi-species *community x invasive snail* experiment below; its true lamp
unit hierarchy is unresolved, so numbers there are descriptive only.
"""
from __future__ import annotations
import argparse, csv, hashlib, itertools, json, math
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA="alan_Xue2026_original6lightingzone_invasive_vs3native_monoculture_v1"
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"database/alan_experiments/xue_sources"
SRC={
 "zone":"8fc7dbdccf92a65ba9e76b2d3d90aebdc836f218",
 "effects":"dcd6f8aa5c027607414523efe4d25e3685f74625",
 "covariance":"773d8b15a31df2770f3bb1758987ca839795c5d8",
 "community":"8ebf5148126ca2610532f1d2cae2a8f51e505165",
 "manifest":"4fcf2d54e8597bb0024fab4bc30b51fb6b694974",
}
ORIGINAL_SHA="8a18f7a30fe8e4656495db511f4e19b0f7e97657233111ba175c7c4007ed4f39"
SPECIES=("Hydrilla_verticillata","Myriophyllum_aquaticum",
         "Myriophyllum_spicatum","Vallisneria_natans")
INV="Myriophyllum_aquaticum"
NATIVE=tuple(x for x in SPECIES if x!=INV)
SEED=20261120
BOOT=5999

def require(ok,msg):
    if not ok:raise RuntimeError("XUE_2026_ORIGINAL_ZONE_SOURCE_STOP: "+msg)

def git_blob(raw):
    return hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()

def read_source():
    frames={}
    names={
       "zone":"monoculture_zone_units.csv",
       "effects":"monoculture_effects.csv",
       "covariance":"monoculture_covariance.csv",
       "community":"community_descriptive_only.csv",
       "manifest":"source_manifest.csv",
    }
    for label,name in names.items():
        data=(BASE/name).read_bytes()
        require(git_blob(data)==SRC[label],"ORIGINAL_PIN_CHANGED: "+name)
        frames[label]=pd.read_csv(BASE/name)
    files=frames["manifest"]
    src=files[files.filename.eq("Monoculture_Growth.xlsx")]
    require(len(src)==1 and src.iloc[0]["sha256"]==ORIGINAL_SHA,
            "Monoculture_Growth.xlsx original independent source SHA changed")
    return frames

def validate_zones(data):
    req={"Light","Zone","Species","biomass_g","n_surviving_cylinders"}
    require(req.issubset(data),"missing original source zone fields")
    require(len(data)==24 and not data.duplicated(["Zone","Species"]).any(),
            "expected 24 unique actual lighting zone x species rows")
    require(set(data.Species)==set(SPECIES),"4 original living macrophyte taxa changed")
    require(data.Zone.nunique()==6,"6 original physical zones changed")
    require(set(data.Light)=={"ALAN","No ALAN"},"true source illumination labels changed")
    counts=data.groupby("Zone").Light.agg(["nunique","first"])
    require(counts["nunique"].eq(1).all(),"original source lighting allocation conflicts")
    require(counts["first"].value_counts().to_dict()=={"ALAN":3,"No ALAN":3},
            "true original 3 lit / 3 dark light-zone allocation changed")
    require((data.n_surviving_cylinders>=2).all(),
            "source missing plants cannot become zero or single-survivor zone")
    require(int(data.n_surviving_cylinders.sum())==77,
            "source original survivor cylinders 77 changed")
    require(np.isfinite(data.biomass_g).all() and (data.biomass_g>0).all(),
            "positive log biomass contrasts require genuine positive source data")
    p=data.pivot(index="Zone",columns="Species",values="biomass_g")
    require(p.shape==(6,4) and p.notna().all().all(),"full original source grid incomplete")
    labels=counts["first"].loc[p.index]
    return p,labels

def index(p):
    return np.log(p[INV].to_numpy(float)) - np.log(p[list(NATIVE)].to_numpy(float)).mean(axis=1)

def evaluate_labels(p,labels,reps=BOOT):
    values=index(p)
    a=np.array(labels.eq("ALAN"))
    require(a.sum()==3 and (~a).sum()==3,"true original treatment count invalid")
    observed=float(values[a].mean()-values[~a].mean())
    combos=list(itertools.combinations(range(len(values)),int(a.sum())))
    require(len(combos)==20,"exact six-unit label allocations not 20")
    allnull=[]
    for ids in combos:
        flag=np.zeros(6,dtype=bool);flag[list(ids)]=True
        allnull.append(float(values[flag].mean()-values[~flag].mean()))
    # Absolute two-sided conditional null p. If original zone assignments
    # aren't source-confirmed randomized, it is a diagnostic, not causal p.
    p_two=sum(abs(v)>=abs(observed)-1e-12 for v in allnull)/len(allnull)
    rng=np.random.default_rng(SEED)
    treated=values[a]
    dark=values[~a]
    draws=np.array([
      rng.choice(treated,len(treated),replace=True).mean()-
      rng.choice(dark,len(dark),replace=True).mean()
      for _ in range(reps)])
    loo=[]
    for i in range(3):
        for j in range(3):
            loo.append(float(np.delete(treated,i).mean()-np.delete(dark,j).mean()))
    return {
      "method":"per-original-lighting-zone log invasive biomass minus log geometric mean of 3 native monoculture biomasses",
      "independent_lit_light_zones":3,
      "independent_dark_light_zones":3,
      "observed_lit_minus_dark_invasive_relative_to_native_log_ratio":observed,
      "invasive_to_native_relative_response_ratio":float(np.exp(observed)),
      "source_lit_perzone_mean_log_index":float(treated.mean()),
      "source_dark_perzone_mean_log_index":float(dark.mean()),
      "true_source_lit_zone_index_range":[float(treated.min()),float(treated.max())],
      "true_source_dark_zone_index_range":[float(dark.min()),float(dark.max())],
      "all_three_lit_zone_indices_exceed_every_dark_zone_index":
          bool(treated.min()>dark.max()),
      "zone_stratified_bootstrap_draws":reps,
      "zone_bootstrap_descriptive_95pct_CI_log_difference":
         [float(x) for x in np.quantile(draws,[.025,.975])],
      "zone_bootstrap_descriptive_95pct_CI_ratio":
         [float(x) for x in np.exp(np.quantile(draws,[.025,.975]))],
      "exact_exchangeable_three_of_six_zone_assignments":len(combos),
      "exact_zone_label_two_sided_p_if_exchangeable":float(p_two),
      "not_proven_original_randomization_p":True,
      "leave_one_lit_and_one_dark_zone_out_log_diff_range":
          [float(min(loo)),float(max(loo))],
      "source_zone_index_by_group_descriptive_only":{
         "ALAN":[float(x) for x in treated],
         "No ALAN":[float(x) for x in dark]},
      "monoculture_not_direct_species_competition":True,
      "not_wet_dry_species_moisture_moderator":True,
    }

def species_effects(p,labels,original_effects,covariance):
    lit=p.loc[labels.eq("ALAN")]
    dark=p.loc[labels.eq("No ALAN")]
    logrr={x:float(np.log(lit[x].mean()/dark[x].mean())) for x in SPECIES}
    for z in original_effects.itertuples(index=False):
        name=str(z.species).replace(" ","_")
        require(name in SPECIES,"unexpected original source effect species")
        require(abs(logrr[name]-float(z.lnRR))<1e-11,
                "stored source original lnRR does not reproduce six zone means")
        require(z.independent_unit=="lighting_zone" and
                z.aggregation=="equal_zone_weight_surviving_cylinder_means",
                "source unit definitions drift")
    w=np.array([(-1/3 if sp!=INV else 1.0) for sp in SPECIES])
    effects=original_effects.set_index("species")
    require(len(effects)==4,"4 original source effects")
    maps={str(x).replace(" ","_"):i for i,x in enumerate(SPECIES)}
    source=np.zeros((4,4),float)
    for z in covariance.itertuples(index=False):
        i=str(z.effect_id_i).replace("xue2026_monoculture_","")
        j=str(z.effect_id_j).replace("xue2026_monoculture_","")
        require(i in maps and j in maps,"covariance species source mismatch")
        source[maps[i],maps[j]]=float(z.covariance)
    require(np.allclose(source,source.T,atol=1e-12) and
            np.linalg.eigvalsh(source).min()>-1e-9,"original source covariance invalid")
    contrast=logrr[INV]-np.mean([logrr[x] for x in NATIVE])
    variance=float(w@source@w)
    require(variance>0,"source covariance invalid")
    return {
      "source_species_monoculture_lit_over_dark_log_mean_ratios":logrr,
      "source_natives_original_log_mean_ratio_average":
        float(np.mean([logrr[x] for x in NATIVE])),
      "source_invasive_lnRR_minus_native_average_lnRR":float(contrast),
      "source_invasive_vs_native_ratio_of_mean_response_ratios":
        float(np.exp(contrast)),
      "within_source_lighting_zone_covariance_delta_approx_variance":variance,
      "within_source_lighting_zone_covariance_delta_approx_SE":float(np.sqrt(variance)),
      "different_from_mean_of_log_zone_indices":True,
      "no_normal_approx_test_claimed_with_only_three_zones_per_arm":True
    }

def community_description(rows):
    require(len(rows)==12 and set(rows.endpoint)=={
        "Native_total_biomass","Invasive_total_biomass"},
        "original 2x3x2 descriptive community source changed")
    require(rows.eligible_for_meta.eq(False).all() and
            rows.independent_lighting_n.astype(str).eq("unresolved").all(),
            "source community real independent fixture status changed")
    def mean(light,herb,endpoint):
        r=rows[(rows.light==light)&(rows.herbivore==herb)&(rows.endpoint==endpoint)]
        require(len(r)==1,"original source community treatment cell missing")
        return float(r.iloc[0].mean_g)
    summary={}
    for l in ("ALAN","No ALAN"):
        original=mean(l,"No snail","Native_total_biomass")
        eaten=mean(l,"Invasive snail","Native_total_biomass")
        invasive=mean(l,"Invasive snail","Invasive_total_biomass")
        invasive_no=mean(l,"No snail","Invasive_total_biomass")
        summary[l]={
          "native_biomass_mean_no_snail_g":original,
          "native_biomass_mean_invasive_snail_g":eaten,
          "native_retention_relative_to_same_light_no_snail":eaten/original,
          "native_suppression_relative_to_same_light_no_snail":1-eaten/original,
          "invasive_biomass_mean_invasive_snail_g":invasive,
          "invasive_biomass_mean_no_snail_g":invasive_no,
          "invasive_retention_relative_to_same_light_no_snail":invasive/invasive_no,
        }
    ratio=summary["ALAN"]["native_retention_relative_to_same_light_no_snail"]/summary["No ALAN"]["native_retention_relative_to_same_light_no_snail"]
    return {
      "summary_by_original_light_arm":summary,
      "descriptive_lit_vs_dark_ratio_of_native_retention_under_invasive_snail":ratio,
      "descriptive_ln_ratio_of_retention_ratios":float(np.log(ratio)),
      "original_source_bucket_count":int(rows.observed_buckets.max()),
      "underlying_real_independent_lamp_count_UNRESOLVED":True,
      "source_animal_grazer_induced_losses_not_true_plant_pollination":True,
      "do_not_infer_statistical_significance_from_independent_bucket_counts":True,
      "original_selection_3_missing_containers_all_ALAN_native_snail":True,
      "not_same_experiment_as_monoculture_50day_zone_outcome":True,
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",required=True)
    args=a.parse_args()
    originals=read_source()
    p,l=validate_zones(originals["zone"])
    ans={
       "schema":SCHEMA,
       "original_study_DOI":"10.1002/ecy.70458",
       "original_ScienceDB_dataset_DOI":"10.57760/sciencedb.35715",
       "original_source_Monoculture_Growth_xlsx_SHA256":ORIGINAL_SHA,
       "original_git_source_pins":SRC,
       "original_3lit_3dark_77surviving_cylinders":True,
       "independent_experimental_light_units":6,
       "same_original_four_species_24_zone_means":True,
       "real_original_4_monoculture_species_names":list(SPECIES),
       "original_invasive_species":INV,
       "original_native_species":list(NATIVE),
       "primary_relative_source_zone_result":evaluate_labels(p,l),
       "secondary_source_covariance_lnRR_result":
           species_effects(p,l,originals["effects"],originals["covariance"]),
       "different_mixed_community_no_source_light_zone_id":
           community_description(originals["community"]),
       "original_source_unobserved_or_dead_containers_not_imputed_zero":True,
       "no_binary_moisture_wet_vs_drier_moderator_in_four_aquatic_species":True,
       "not_many_independent_cross_study_plant_light_experiments":True,
       "no_direct_seed_or_pollination_measured":True,
       "not_a_verified_mechanism_for_discovery_era_macro_moist_decliners":True,
    }
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(ans,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    m=ans["primary_relative_source_zone_result"]
    n=ans["secondary_source_covariance_lnRR_result"]
    c=ans["different_mixed_community_no_source_light_zone_id"]
    print("XUE2026_ORIGINAL_ZONE_WISE_INVADER_NATIVE_RESPONSE "+
          json.dumps({
            "zone_index":m,
            "original_species_lit_dark":n,
            "original_grazed_community_descriptive_only":c
          },sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
