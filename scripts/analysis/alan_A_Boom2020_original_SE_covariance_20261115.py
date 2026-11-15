#!/usr/bin/env python3
"""Boom 2020 original JATS Table 1: why per-arm SE cannot yield a paired meta SE.

Author-reported total *estimated* seeds/flower (formula) ± source arm SE:
dark117±40, white86±28, green61±20, red90±28, four shared original sites.
Re-fetch DOI-verified Europe PMC primary JATS; preserve source point and SE.
Use Cauchy-Schwarz covariance bounds, NOT a new fitted CI or pairwise p.
"""
from __future__ import annotations
import hashlib,json,math,re,urllib.request
from pathlib import Path
import xml.etree.ElementTree as ET
from scripts.analysis.alan_original_Boom2020_nursery_pollinator_seedpredation_spectra_20261114 import (
    XML_URL,ARTICLE_DOI,SOURCE_NAMES,FIELD_SITES,get_table,normalized_text,decomposition
)
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"results/alan_A_Boom2020_estimated_seed_output_pairwise_SE_covariance_bounds_20261115.json"
EXPECTED_SE={"Dark":40.0,"White":28.0,"Green":20.0,"Red":28.0}
EXPECTED_MEAN={"Dark":117.0,"White":86.0,"Green":61.0,"Red":90.0}
def exact_source():
    req=urllib.request.Request(XML_URL,headers={"User-Agent":"ALAN-Boom2020-paired-covariance-source/20261115"})
    with urllib.request.urlopen(req,timeout=60) as r:blob=r.read(4_000_000)
    root=ET.fromstring(blob)
    table=get_table(root)
    wrap=None
    for x in root.iter("table-wrap"):
        cap=normalized_text(x.find("caption")) if x.find("caption") is not None else ""
        if "Proportion" in cap and "ovules" in cap and ("infested" in cap or "fruits" in cap):
            if wrap is not None:raise ValueError("Multiple original Table 1 candidates")
            wrap=x
    if wrap is None:raise ValueError("No author Table 1")
    printed={}
    for tr in wrap.iter("tr"):
        cells=[normalized_text(td) for td in tr if td.tag.rsplit("}",1)[-1] in ("td","th")]
        if not cells or cells[0] not in SOURCE_NAMES:continue
        # Preserve original non-pairwise per-arm reported standard errors.
        m=re.search(r"(?P<mean>[0-9]+)\s*(?:±|\+/-|\+-)\s*(?P<se>[0-9]+)",cells[4])
        if m is None:raise ValueError("Original estimated-seed mean±SE missing in "+cells[0])
        printed[cells[0]]={"estimated_seed_mean":float(m.group("mean")),
                           "source_arm_SE":float(m.group("se")),
                           "original_source_estimated_not_observed_seed_count":True}
        if printed[cells[0]]["estimated_seed_mean"]!=EXPECTED_MEAN[cells[0]]:
            raise ValueError("Original seeds/flower source mean differs")
        if printed[cells[0]]["source_arm_SE"]!=EXPECTED_SE[cells[0]]:
            raise ValueError("Original seeds/flower source SE differs")
        if table[cells[0]]["estimated_seeds_per_flower_printed"]!=EXPECTED_MEAN[cells[0]]:
            raise ValueError("Original source Table 1 stage reconciliation failed")
    if set(printed)!=set(SOURCE_NAMES):raise ValueError("Missing source arm")
    return printed,table,hashlib.sha256(blob).hexdigest()

def covariance_bounds(mean_dark,se_dark,mean_treat,se_treat):
    """Var(treat-dark)=s_t^2+s_d^2 -2 rho*s_t*s_d, -1<=rho<=1.

    These are logical bounds from per-arm SEs; not a valid confidence interval.
    """
    if se_dark<0 or se_treat<0 or mean_dark<=0 or mean_treat<=0:
        raise ValueError("Source mean/SE invalid")
    delta=mean_treat-mean_dark
    min_se=abs(se_dark-se_treat)
    hypothetical_max_abs_t_4site=(abs(delta)/min_se if min_se>0 else None)
    return {
        "reported_mean_treatment_minus_dark":mean_treat-mean_dark,
        "hypothetical_best_case_abs_t_if_four_independent_matched_sites":hypothetical_max_abs_t_4site,
        "hypothetical_4site_df3_critical_t_95_two_sided":3.182446,
        "hypothetical_best_case_t_exceeds_df3_threshold":
             hypothetical_max_abs_t_4site>3.182446 if hypothetical_max_abs_t_4site is not None else None,
        "site_t_computation_is_not_original_paper_test":True,
        "reported_mean_treatment_to_dark":mean_treat/mean_dark,
        "possible_pair_difference_SE_min_if_rho_plus1":abs(se_dark-se_treat),
        "possible_pair_difference_SE_max_if_rho_minus1":se_dark+se_treat,
        "illustrative_independent_arm_SE_if_rho_zero":math.hypot(se_dark,se_treat),
        "published_covariance_rho_not_available":True,
        "not_a_computed_meta_lnRR_sampling_variance":True,
        "original_full_4_treatment_ANOVA_p":0.63
    }

def calculate(printed,table,jats_sha):
    if set(printed)!=set(SOURCE_NAMES):raise AssertionError("Four published arms required")
    if set(table)!=set(SOURCE_NAMES):raise AssertionError("Wrong source")
    dark=printed["Dark"];contrast={}
    for light in ("White","Green","Red"):
        v=printed[light]
        contrast[light]=covariance_bounds(dark["estimated_seed_mean"],
            dark["source_arm_SE"],v["estimated_seed_mean"],v["source_arm_SE"])
        ratio=decomposition(table)[light]["ratio_to_dark_recomputed"]
        contrast[light]["ratio_from_rounded_fruit_infestation_ovule_stage_means"]=ratio
        contrast[light]["stage_product_and_published_seed_mean_differ_by_rounding"]=True
        contrast[light]["source_shares_dark_plant_sites"]=True
    return {
      "schema":"ALAN_A_Boom2020_original_estimated_seeds_per_flower_covariance_identifiability_v1",
      "doi":ARTICLE_DOI,"source_JATS":XML_URL,"source_jats_sha256":jats_sha,
      "source_field_sites":FIELD_SITES,"original_source_arms":printed,
      "source_4arm_estimated_seed_output_ANOVA_p":0.63,
      "contrasts_vs_shared_dark":contrast,
      "limits":[
         "Paper's per-arm seed output is a formula product, not a measured original-flower individual seed count",
         "Published ±SE per arm does not provide covariance of repeated site-level treatments",
         "Four-field site paired covariance unknown; SE bounds are arithmetic, not confidence intervals or estimates",
         "All three spectrum contrasts reuse the same original dark control: one paper, not three experiments",
         "Original author 4-treatment estimated-seed-output ANOVA p=0.63, not a significant global effect",
         "Seed/ovule fertilization measured conditional on uninfested fruit, selected after pollinator-larval damage",
         "No common light/dark open total seed lnRR cluster variance or independent dark-night V recovered",
         "Never convert individual 600+ flower counts to independent 600+ light installations"
      ],
      "global_ALAN_night_dependence_beta_V":None
    }

def main():
    printed,table,sha=exact_source()
    a=calculate(printed,table,sha)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(a,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("ALAN_BOOM_SOURCE_SE_COVARIANCE="+json.dumps(a,ensure_ascii=False,allow_nan=False))
if __name__=="__main__":main()
