#!/usr/bin/env python3
"""Source-original pollinating seed-predator antagonism under ALAN spectra.

Boom et al. 2020, Scientific Reports, DOI 10.1038/s41598-020-75471-1.
Read the publisher's OA JATS XML Table 1 through Europe PMC at CI time, rather
than inventing a local seed-level experiment. Output only treatment-level
published aggregates and an exact decomposition of the ORIGINAL defined
estimated seed-output *formula*. It is NOT causal mediation, not a global
pollination syndrome ALAN effect, and not a test of predator substitution.
"""
from __future__ import annotations
import argparse, hashlib, json, math, re
from datetime import datetime, timezone
from pathlib import Path
import urllib.request
import xml.etree.ElementTree as ET

ARTICLE_DOI="10.1038/s41598-020-75471-1"
PMCID="PMC7591485"
XML_URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7591485/fullTextXML"
SOURCE_NAMES=("Dark","White","Green","Red")
EXPECTED={
    "Dark":  (0.73,0.28,0.56,117,689,497,28348),
    "White": (0.78,0.44,0.49,86,680,518,23822),
    "Green": (0.75,0.53,0.43,61,596,429,29127),
    "Red":   (0.71,0.33,0.47,90,663,514,27650),
}
FIELD_SITES=4
AVERAGE_OVULES=404
REPORTED_CHI2_FRUIT_P=0.17
REPORTED_CHI2_INFEST_P_LESS_THAN=0.001
REPORTED_OVULE_FERTILIZATION_P_LESS_THAN=0.01
REPORTED_ESTIMATED_OUTPUT_ANOVA_P=0.63

def require(ok,msg):
    if not bool(ok):
        raise ValueError("ALAN_BOOM_PRIMARY: "+str(msg))

def normalized_text(el):
    return re.sub(r"\s+"," ","".join(el.itertext())).strip()

def get_table(root):
    """Do not silently substitute a numerical table from another paper."""
    article_id=[normalized_text(x).lower() for x in root.iter("article-id")
                if x.attrib.get("pub-id-type")=="doi"]
    require(ARTICLE_DOI.lower() in article_id,
            "EuropePMC returned wrong DOI instead of Boom 2020")
    tables=[]
    for x in root.iter("table-wrap"):
        a=normalized_text(x.find("caption")) if x.find("caption") is not None else ""
        if "Proportion" in a and "ovules" in a and ("infested" in a or "fruits" in a):
            tables.append(x)
    require(len(tables)==1,"Expected one actual Table 1 showing fruits, infestation and ovules")
    r={}
    for tr in tables[0].iter("tr"):
        cells=[normalized_text(td) for td in tr if td.tag.rsplit("}",1)[-1] in ("td","th")]
        if not cells or cells[0] not in SOURCE_NAMES:
            continue
        require(len(cells)>=5,"Original source row missing four outcome columns")
        treatment=cells[0]
        parts=[]
        for i,col in enumerate(cells[1:4]):
            match=re.search(r"(0(?:\.\d+)?|1(?:\.0+)?)",col)
            nmatch=re.search(r"\(\s*([0-9,]+)\s*\)",col)
            require(bool(match) and bool(nmatch),
                    "No published treatment mean and denominator in "+treatment+" column "+str(i))
            parts.append((float(match.group(1)),int(nmatch.group(1).replace(",",""))))
        src_out=re.search(r"\b([0-9]+)\s*(?:±|\+/-|\+-)",cells[4])
        require(src_out is not None,"Original published estimated seed output absent")
        r[treatment]={
            "fruit_set":parts[0][0],
            "fruit_set_source_nested_flowers":parts[0][1],
            "infestation_given_fruit":parts[1][0],
            "infestation_source_nested_fruits":parts[1][1],
            "ovule_fertilization_given_uninfested_fruit":parts[2][0],
            "ovule_sample_nested_in_fruit_plant_transect_site":parts[2][1],
            "estimated_seeds_per_flower_printed":int(src_out.group(1)),
        }
    require(set(r)==set(SOURCE_NAMES),"Missing verified original 4 light-spectrum arms")
    for name,entry in r.items():
        ref=EXPECTED[name]
        got=(entry["fruit_set"],entry["infestation_given_fruit"],
             entry["ovule_fertilization_given_uninfested_fruit"],
             entry["estimated_seeds_per_flower_printed"],
             entry["fruit_set_source_nested_flowers"],
             entry["infestation_source_nested_fruits"],
             entry["ovule_sample_nested_in_fruit_plant_transect_site"])
        require(got==ref,"Original published Table 1 values changed: "+name)
    return r

def decomposition(data):
    """Exact identity on published rounded treatment means, NOT causal effect.
       E_t = p(fruit)_t * (1-p(seed predator infestation)|fruit)_t
             * p(ovule fertilized|noninfested fruit)_t * 404 ovules/fruit.
    """
    dark=data["Dark"]
    p0=dark["fruit_set"]
    a0=1-dark["infestation_given_fruit"]
    o0=dark["ovule_fertilization_given_uninfested_fruit"]
    baseline=p0*a0*o0*AVERAGE_OVULES
    out={}
    for name in SOURCE_NAMES:
        entry=data[name]
        p=entry["fruit_set"];a=1-entry["infestation_given_fruit"]
        o=entry["ovule_fertilization_given_uninfested_fruit"]
        estimate=p*a*o*AVERAGE_OVULES
        require(abs(estimate-entry["estimated_seeds_per_flower_printed"])<=3.5,
                "Calculated output does not match original printed estimates within rounding")
        factors={"fertilized_flower_factor":p/p0,
                 "escape_seed_predation_factor":a/a0,
                 "fertilized_ovule_factor":o/o0}
        product=math.prod(factors.values())
        require(abs(estimate/baseline-product)<1e-10,"Product identity drift")
        out[name]={
          "published_original_fruit_set":p,
          "published_original_fruit_infestation":entry["infestation_given_fruit"],
          "published_original_ovule_fertilization_in_selected_uninfested_fruits":o,
          "original_nested_denominators":{
             "fruit_flowers":entry["fruit_set_source_nested_flowers"],
             "infested_fruits":entry["infestation_source_nested_fruits"],
             "sampled_ovules_not_independent":entry["ovule_sample_nested_in_fruit_plant_transect_site"]},
          "published_original_estimated_seeds_per_flower":entry["estimated_seeds_per_flower_printed"],
          "recomputed_from_rounded_source_means":estimate,
          "ratio_to_dark_recomputed":estimate/baseline,
          "percent_change_to_dark_from_rounded_means":100*(estimate/baseline-1),
          "algebraic_source_product_ratio_factors":factors,
          "algebraic_source_log_contribution":{
              key:math.log(v) for key,v in factors.items()},
          "not_causal_mediator_effects":True,
        }
    return out

def main(out):
    out=Path(out)
    require(not out.exists(),"Never overwrite previously audited source output")
    req=urllib.request.Request(XML_URL,headers={
        "User-Agent":"ALAN-open-science-nursery-pollination-source/1.0",
        "Accept":"application/xml,text/xml,*/*"})
    with urllib.request.urlopen(req,timeout=35) as resp:
        blob=resp.read(5_000_000)
    require(blob.startswith(b"<?xml") or b"<article" in blob[:300],
            "Original full text source did not return JATS XML")
    root=ET.fromstring(blob)
    data=get_table(root)
    derived=decomposition(data)
    result={
     "schema":"ALAN_original_silene_hadena_mutualism_seedpredation_spectral_20261114_v1",
     "status":"PRIMARY_PUBLISHED_TABLE_VERIFIED_ALGEBRAIC_MECHANISM_ONLY",
     "source":{
       "paper_DOI":ARTICLE_DOI,"PMCID":PMCID,"original_table":"Table 1",
       "primary_fulltext_JATS_url":XML_URL,
       "downloaded_raw_JATS_sha256":hashlib.sha256(blob).hexdigest(),
       "source_year":2020,"independent_field_sites":FIELD_SITES,
       "light_spectra_four_arm":["Dark","White","Green","Red"],
       "mean_ground_lux_original_reported":7.6,
       "mean_ground_lux_source_SE":1.2,
       "mean_ovules_per_fruit_reported":AVERAGE_OVULES,
       "estimated_seed_output_is_formula_not_counted_individual_seeds":True,
       "taxon_of_dual_role_moth":"Hadena bicruris",
       "plant":"Silene latifolia",
       "offspring_develop_in_fruit_and_eat_seeds":True,
       "contemporaneous_nursery_pollination_and_seed_predation_measured_as_same_system":True
     },
     "source_original_reported_inferential_tests":{
       "fruit_set_light_fourarm_Chi2_p":REPORTED_CHI2_FRUIT_P,
       "infestation_light_fourarm_Chi2_p_less_than":REPORTED_CHI2_INFEST_P_LESS_THAN,
       "ovule_fertilized_Chi2_p_less_than":REPORTED_OVULE_FERTILIZATION_P_LESS_THAN,
       "estimated_seeds_per_flower_ANOVA_p":REPORTED_ESTIMATED_OUTPUT_ANOVA_P,
       "do_not_claim_mechanistically_derived_seed_decline_is_significant":True
     },
     "original_light_treatment_table_and_decomposition":derived,
     "functional_interpretation":{
       "pollinator_same_animal_taxon_as_larval_seed_predator":True,
       "ALAN_changes_infested_fruit_rate_without_resolved_fruit_set_change":True,
       "comparison_does_not_isolate_male_vs_female_Hadena_pollination":True,
       "no_independent_actual_other_pollinator_pollen_deposition":True,
       "ovule_success_conditioned_on_uninfested_fruits_may_be_selection_biased":True,
       "original_flowers_fruits_ovules_nested_not_independent_sample_sizes":True,
       "no_actual_recruitment_population_change_or_island_response":True,
       "no_moth_predator_of_adult_floral_visitors_measured":True,
       "not_a_predator_substituting_for_moth_pollen_service":True,
       "no_global_ALAN_by_syndrome_coefficient_estimated":True,
       "spectral_light_geometry_and_local_lux_not_satellite_PCNL":True
     }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_ORIGINAL_NURSERY_POLLINATOR_SEED_PREDATION_SPECTRAL",
          json.dumps({"source":result["source"],
                      "reported_inference":result["source_original_reported_inferential_tests"],
                      "calculated":derived},ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    main(args.out)
