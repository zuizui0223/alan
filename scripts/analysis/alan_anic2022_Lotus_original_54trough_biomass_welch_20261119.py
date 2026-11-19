#!/usr/bin/env python3
"""Source-original 2016 Lotus pedunculatus biomass 3×3 × six physical troughs.

Extracts author Figshare raw numeric PDF (665 species×physical-trough rows)
via SHA-verified prior script; uses only per-cell original means/SD and
physical unit n (not raw rows in outputs). Computes equal-trophic-stratum
raw biomass contrasts and Welch–Satterthwaite approximations. These are
exploratory and unadjusted, NOT authors' GLM/Tukey results and not a
causal plant fitness / food-web mediation model.

Crucial: source PDF light codes L1/L2/L3 and trophic codes T1/T2/T3
are UNMAPPED to unlit/amber/white or plants/herbivores/predators until
verified with an original author treatment key. No guessed mapping.
"""
import argparse
import json
import math
from pathlib import Path
from scipy.stats import t as student_t
from scipy.stats import f as fisher_f

from scripts.analysis.alan_anic2022_original_pdf_biomass_54trough_parse_20261119 import (
    FILE_SHA256,FILE_ID,FIGSHARE_API,get,read_original_biomass_pdf
)

TARGET="Lotus pedunculatus"
LIGHT=("L1","L2","L3")
TROPHIC=("T1","T2","T3")
SCHEMA="alan_20261119_Anic_original_lotus_54trough_equal_trophic_light_contrasts_v1"


def layer(rows):
    raw={(x["original_light_code"],x["original_trophic_code"]):x for x in rows
         if x["species"]==TARGET}
    if set(raw)!={(l,t) for l in LIGHT for t in TROPHIC}:
        raise ValueError("Primary original Lotus species does not survive in every randomization cell")
    out={}
    for l in LIGHT:
        cells=[raw[(l,t)] for t in TROPHIC]
        if any(c["original_mesocosms_with_this_species"]!=6 for c in cells):
            raise ValueError("Original physically independent Lotus trough n must be six per cell")
        out[l]={
            "original_source_light_code":l,
            "source_physical_mesocosms":sum(c["original_mesocosms_with_this_species"] for c in cells),
            "equal_trophic_cell_mean_biomass_g":sum(c["original_biomass_total_g_arithmetic_mean"] for c in cells)/3,
            "trophic_cell_level_descriptives":[
                {"trophic_code":t,
                 "n_original_physical_mesocosms":raw[(l,t)]["original_mesocosms_with_this_species"],
                 "mean_biomass_g":raw[(l,t)]["original_biomass_total_g_arithmetic_mean"],
                 "sample_SD_biomass_g":raw[(l,t)]["original_biomass_total_g_sample_SD"]}
                for t in TROPHIC
            ],
        }
    return raw,out


def balanced_original_factorial_anova(raw):
    """Descriptive classical balanced original light×trophic ANOVA.

    Using true original cell n, sample mean and sample SD is algebraically
    equivalent to a 3×3 OLS two-way ANOVA on the 54 original physical troughs
    under the *homoskedastic* residual assumption. This is a DIFFERENT
    fixed-model estimand from authors' AIC selected GLMs and Tukey contrasts.
    """
    if len(raw)!=9 or any(c["original_mesocosms_with_this_species"]!=6 for c in raw.values()):
        raise ValueError("Original 3×3×6 balanced factorial source required")
    n=6
    cellmean={(l,t):raw[(l,t)]["original_biomass_total_g_arithmetic_mean"]
              for l in LIGHT for t in TROPHIC}
    mean_l={l:sum(cellmean[(l,t)] for t in TROPHIC)/3 for l in LIGHT}
    mean_t={t:sum(cellmean[(l,t)] for l in LIGHT)/3 for t in TROPHIC}
    grand=sum(cellmean.values())/9
    ss_light=n*len(TROPHIC)*sum((mean_l[l]-grand)**2 for l in LIGHT)
    ss_trophic=n*len(LIGHT)*sum((mean_t[t]-grand)**2 for t in TROPHIC)
    ss_interaction=n*sum(
        (cellmean[(l,t)]-mean_l[l]-mean_t[t]+grand)**2
        for l in LIGHT for t in TROPHIC)
    ss_error=sum(
        (n-1)*raw[(l,t)]["original_biomass_total_g_sample_SD"]**2
        for l in LIGHT for t in TROPHIC)
    df_error=len(LIGHT)*len(TROPHIC)*(n-1)
    mse=ss_error/df_error
    if mse<=0:raise ValueError("Original biomass residual MSE invalid")
    effects={}
    for key,ss,df in [
        ("LIGHT",ss_light,len(LIGHT)-1),
        ("TROPHIC",ss_trophic,len(TROPHIC)-1),
        ("LIGHT_X_TROPHIC",ss_interaction,(len(LIGHT)-1)*(len(TROPHIC)-1))
    ]:
        F=(ss/df)/mse
        effects[key]={"SS_g2":ss,"df":df,"F":F,
                      "unadjusted_parametric_p_homoskedastic":fisher_f.sf(F,df,df_error)}
    return {
        "design":"unadjusted balanced fixed 3×3 factorial ANOVA on published original physical-trough means and within-cell sample SD",
        "independently_assigned_physical_trough_n":n*len(raw),
        "physical_troughs_per_cell":n,
        "error_df":df_error,
        "source_residual_within_cell_MSE_g2":mse,
        "primary_original_species":TARGET,
        "factors":effects,
        "not_author_original_AIC_selected_GLM":True,
        "homoskedastic_assumption_unverified":True,
        "not_corrected_across_13_species":True,
        "no_independent_2018_herbivore_source_joint_fit":True,
        "not_a_reproductive_seed_fitness_effect":True,
    }


def source_contrasts(d,primary=TARGET):
    if d["source_original_figshare_sha256"]!=FILE_SHA256:
        raise ValueError("Original author Figure numeric PDF source SHA drift")
    if not d["experimental_light_trophic_full_factorial_verified"]:
        raise ValueError("Cannot fit numeric study contrast without validated 54-trough factorial")
    if d["original_unique_physical_mesocosms_in_biomass_table"]!=54:
        raise ValueError("Original 54 independently assigned experimental units not recovered")
    if d["original_species_rows_parsed"]!=665:
        raise ValueError("Original 665 observed species×trough records drifted")
    raw,out=layer(d["source_original_species_by_light_trophic_cell_descriptive_means"])
    comparisons=[]
    for l0,l1 in [("L1","L3"),("L1","L2"),("L2","L3")]:
        within=[]
        for trophic in TROPHIC:
            a,b=raw[(l0,trophic)],raw[(l1,trophic)]
            n0=a["original_mesocosms_with_this_species"]
            n1=b["original_mesocosms_with_this_species"]
            m0,m1=a["original_biomass_total_g_arithmetic_mean"],b["original_biomass_total_g_arithmetic_mean"]
            s0,s1=a["original_biomass_total_g_sample_SD"],b["original_biomass_total_g_sample_SD"]
            if min(n0,n1)<2 or min(s0,s1)<0:
                raise ValueError("No physical replicate variance for one original experimental treatment")
            v0=s0*s0/n0;v1=s1*s1/n1
            v=v0+v1
            df=v*v/(v0*v0/(n0-1)+v1*v1/(n1-1)) if v>0 else None
            within.append({
                "original_trophic_code":trophic,
                "original_physical_mesocosms_per_light_arm":[n0,n1],
                "raw_treatment_means_g":[m0,m1],
                "source_unadjusted_mean_difference_g":m1-m0,
                "sampling_variance_of_difference_g2":v,
                "Welch_small_sample_df":df,
            })
        # Marginal light difference averaged equally across the
        # 3 randomized trophic strata: 6 physical troughs per cell.
        diff=sum(c["source_unadjusted_mean_difference_g"] for c in within)/3
        v=sum(c["sampling_variance_of_difference_g2"] for c in within)/9
        se=math.sqrt(v)
        df=(sum(c["sampling_variance_of_difference_g2"] for c in within)**2/
            sum(c["sampling_variance_of_difference_g2"]**2/c["Welch_small_sample_df"]
                for c in within if c["Welch_small_sample_df"]))
        critical=student_t.ppf(.975,df)
        p=2*student_t.sf(abs(diff/se),df)
        ref=out[l0]["equal_trophic_cell_mean_biomass_g"]
        cmp=out[l1]["equal_trophic_cell_mean_biomass_g"]
        comparisons.append({
            "original_unmapped_light_code_contrast":f"{l1} minus {l0}",
            "reference_light_code":l0,"comparator_light_code":l1,
            "original_paper_species":primary,
            "source_original_trophic_strata":3,
            "original_physical_mesocosms_each_light_arm":18,
            "reference_original_biomass_mean_g":ref,
            "comparator_original_biomass_mean_g":cmp,
            "comparator_over_reference_unadjusted_arithmetic_ratio":cmp/ref if ref>0 else None,
            "difference_g":diff,
            "independent_trough_Welch_approx_SE_g":se,
            "Welch_Satterthwaite_df":df,
            "exploratory_unadjusted_t_P":p,
            "exploratory_unadjusted_95pct_CI_g":[diff-critical*se,diff+critical*se],
            "original_trophic_condition_cell_contrasts":within,
            "model":"equal_trophic_strata independent-trough raw mean difference and Welch-Satterthwaite t; no covariate adjustment or multi-comparison correction",
            "not_original_author_AIC_GLM_or_Tukey_p":True,
            "source_original_light_codes_not_mapped_to_spectral_treatments":True,
            "not_cross_study_meta_analytical_effect":True,
        })
    order=sorted(range(len(comparisons)),key=lambda i:comparisons[i]["exploratory_unadjusted_t_P"])
    adjusted=[0]*len(comparisons)
    mx=0
    for rank,idx in enumerate(order):
        rawp=comparisons[idx]["exploratory_unadjusted_t_P"]
        mx=max(mx,min(1.,rawp*(len(comparisons)-rank)))
        adjusted[idx]=mx
    for c,adj in zip(comparisons,adjusted):
        c["exploratory_Holm_multiplicity_adjusted_P_three_light_pairs"]=adj
    factorial=balanced_original_factorial_anova(raw)
    return {
        "schema":SCHEMA,
        "original_2022_primary_DOI":"10.1002/ece3.9157",
        "original_data_Figshare_PDF_SHA256":FILE_SHA256,
        "original_source_species_observed_total":13,
        "original_species_by_trough_actual_raw_records":665,
        "original_independent_physical_light_trophic_troughs":54,
        "unique_original_programmes":1,
        "primary_species":primary,
        "selected_due_to_original_published_species_specific_biomass_response":True,
        "original_source_codes_preserved_without_spectral_identity":True,
        "original_light_level_aggregate_physical_trough_descriptives":out,
        "source_original_exploratory_light_pair_contrasts":comparisons,
        "original_54trough_Lotus_balanced_3x3_factorial_ANOVA":factorial,
        "source_original_light_code_mapping_verified_to_unlit_amber_white":False,
        "source_2018_Dryad_herbivore_rows_joined_by_original_trough":False,
        "no_reproductive_seed_or_fruit_endpoint_measured":True,
        "exploratory_analytic_multiple_testing_and_original_author_pvalue_separate":True,
        "no_original_plot_IDs_or_species_by_trough_values_exported":True,
        "source_original_author_GLM_reanalysis_not_fitted":True,
        "no_predator_mediation_fitted":True,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    metadata=json.loads(get(FIGSHARE_API))
    files=[x for x in metadata.get("files",[]) if x.get("id")==FILE_ID]
    if len(files)!=1:
        raise ValueError("Original source Figshare file changed")
    b=get(files[0]["download_url"])
    verified=read_original_biomass_pdf(b,strict=True)
    ans=source_contrasts(verified)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(ans,indent=2,ensure_ascii=False)+"\n")
    for c in ans["source_original_exploratory_light_pair_contrasts"]:
        print("ANIC2022_SOURCE_LOTUS_ORIGINAL_LIGHT_CODE_CONTRAST",
              c["original_unmapped_light_code_contrast"],
              "delta_g",round(c["difference_g"],4),
              "SE",round(c["independent_trough_Welch_approx_SE_g"],4),
              "CI",c["exploratory_unadjusted_95pct_CI_g"],
              "P_nominal",c["exploratory_unadjusted_t_P"],
              "P_Holm_3",c["exploratory_Holm_multiplicity_adjusted_P_three_light_pairs"])
    print("ANIC2022_LOTUS_ORIGINAL_54TROUGH_FACTORIAL_ANOVA",ans["original_54trough_Lotus_balanced_3x3_factorial_ANOVA"])
    print("ANIC2022_LOTUS_TRUE_SOURCE_54_TROUGH_CONTRAST_GATE",
          ans["source_original_light_code_mapping_verified_to_unlit_amber_white"])


if __name__=="__main__":
    main()
