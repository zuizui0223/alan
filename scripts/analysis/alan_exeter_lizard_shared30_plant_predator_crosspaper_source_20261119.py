#!/usr/bin/env python3
"""Cross-publication physical-experiment identity and original archive access.

The published 30 physical 4×4m Lizard grassland light installations are shared
by Bennie2018 plant phenology and Davies2017 spider/beetle sampling. The two
papers are 2 endpoints, NOT 2 independently randomized lighting experiments.

Source archive scans identify public raw-file availability only, with no
plot-level data export or fictitious joins by treatment arm alone.
"""
import argparse
from collections import defaultdict
from html.parser import HTMLParser
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

DEFAULT=Path("results/all_species_forest/biotic_interactions/alan_exeter_lizard_shared_30lightplots_crosspaper_provenance_20261119.json")
SCHEMA="alan_exeter_lizard_2012_shared_30_light_plot_programme_cross_trophic_20261119_v1"
MAX_BYTES=1_200_000

class LinkReader(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links=[]
    def handle_starttag(self, tag, attrs):
        if tag!="a":return
        href=dict(attrs).get("href")
        if isinstance(href,str):self.links.append(href)

def inspect_html_links(html,base):
    q=LinkReader()
    q.feed(html.decode("utf-8","replace"))
    found=[]
    for src in q.links:
        resolved=urllib.parse.urljoin(base,src)
        url=urllib.parse.urlsplit(resolved)
        if url.scheme not in ("http","https"):
            continue
        if url.hostname not in ("ore.exeter.ac.uk","repository.exeter.ac.uk"):
            continue
        if ("/bitstream/" in url.path or "/download" in url.path or
            url.path.lower().endswith((".csv",".tsv",".xlsx",".xls",".rdata",".rds",".zip"))):
            if resolved not in found:found.append(resolved)
    return found[:40]

def read_bounded(url):
    req=urllib.request.Request(url,headers={"User-Agent":"ALAN-source-authenticity-Exeter-2017-2018/1.0",
                                           "Accept":"text/html,application/json,*/*"})
    with urllib.request.urlopen(req,timeout=25) as response:
        body=response.read(MAX_BYTES+1)
        if len(body)>MAX_BYTES:raise ValueError("Source response over size limit")
        return body,response.geturl(),response.headers.get("Content-Type","")

def source_scan(links=inspect_html_links, reader=read_bounded):
    endpoints={
      "Bennie2018_flowering_original_archived_handle":[
        "https://hdl.handle.net/10871/27138",
        "https://ore.exeter.ac.uk/repository/handle/10871/27138",
        "https://ore.exeter.ac.uk/repository/server/api/core/handles/10871/27138"
      ],
      "Davies2017_spider_beetle_article_archived_handle":[
        "https://hdl.handle.net/10871/25044",
        "https://ore.exeter.ac.uk/repository/handle/10871/25044",
        "https://ore.exeter.ac.uk/repository/server/api/core/handles/10871/25044"
      ]
    }
    result={}
    for label,urls in endpoints.items():
        checked=[]
        for url in urls:
            try:
                body,final,typ=reader(url)
                finalhost=urllib.parse.urlsplit(final).hostname or ""
                entry={"status":"SOURCE_METADATA_RETRIEVED",
                       "requested_url":url,
                       "final_source_host":finalhost,
                       "content_type":typ,
                       "metadata_bytes":len(body)}
                if "html" in typ.lower() or body.lstrip().startswith(b"<!DOCTYPE"):
                    entry["official_archived_file_candidates"]=links(body,final)
                elif "json" in typ.lower() or body.lstrip().startswith(b"{"):
                    # Preserve status and metadata size, do not infer treatment
                    # plot IDs or source field data from item metadata.
                    entry["official_archived_file_candidates"]=[]
                checked.append(entry)
            except (ValueError,urllib.error.URLError,urllib.error.HTTPError,TimeoutError) as e:
                checked.append({
                    "status":"SOURCE_STOP","requested_url":url,
                    "error_class":type(e).__name__,
                    **({"http_status":e.code} if isinstance(e,urllib.error.HTTPError) else {})
                })
        result[label]=checked
    return {
        "schema":"alan_20261119_exeter_two_archives_original_30plot_joint_source_access_v1",
        "field_experimental_programme":"ExeterLizardLED30_2012_2015",
        "physical_independently_assigned_light_plots":30,
        "original_lighting_treatment_arms":5,
        "joint_original_plant_arachnid_plot_records_verified":False,
        "verified_source_joint_multi_trophic_treatment_effect_estimated":False,
        "public_archival_page_access":result,
        "no_original_plot_ID_or_arthropod_row_exported":True,
        "metadata_success_not_raw_multitrophic_join":True,
        "same_field_2015_exeter_mesocosm54_not_part_of_this_programme":True
    }


def validate(registry):
    if registry.get("schema")!=SCHEMA or registry.get("programme_id")!="ExeterLizardLED30_2012_2015":
        raise ValueError("Wrong original published light installation programme")
    site=registry["site"]
    if site["original_independently_randomized_light_plots"]!=30:
        raise ValueError("Original physical light-unit count drifted")
    if site["physical_plot_area_m2"]!=16 or site["shared_field_years"]!=[2012,2013,2014]:
        raise ValueError("Original 4×4m 2012–14 overlap not source-matched")
    arms=registry["physical_treatment_arms"]
    if len(arms)!=5 or set(x["canonical"] for x in arms)!={"CONTROL","HIW","LIW","LIWT","AMB"}:
        raise ValueError("Physical five-light treatment arms missing")
    if sum(x["replicate_physical_plots"] for x in arms)!=30:
        raise ValueError("Source treatment arms are not six independent plots each")
    for x in arms:
        if x["replicate_physical_plots"]!=6:
            raise ValueError("Mismatched source treatment light plot repetitions")
        if not x.get("Bennie2018") or not x.get("Davies2017"):
            raise ValueError("Source same-arm alias missing from one paper")
    expected={"CONTROL":"unlit control","HIW":"HIW","LIW":"DW","LIWT":"DWT","AMB":"AMB"}
    if {x["canonical"]:x["Davies2017"] for x in arms}!=expected:
        raise ValueError("Original cross-paper source lighting aliases require re-audit")
    pubs=registry["publications"]
    if len(pubs)!=2 or set(x["doi"] for x in pubs)!={"10.1111/1365-2664.12927","10.1111/gcb.13615"}:
        raise ValueError("Original twin scientific paper sources changed")
    bydoi={x["doi"]:x for x in pubs}
    a=bydoi["10.1111/1365-2664.12927"]
    b=bydoi["10.1111/gcb.13615"]
    if not set(site["shared_field_years"]).issubset(set(a["original_years"]) & set(b["original_years"])):
        raise ValueError("Cannot join unrelated site-years")
    if len(a["original_species_analyzed_flowering"])!=3:
        raise ValueError("Only original three grasses had modeled flowering dates")
    if b["source_arthropods_caught"]!=5180 or b["source_arthropod_taxa"]!=136:
        raise ValueError("Published original ground predator source stats changed")
    e=registry["source_linkage_evidence"]
    if not e["physical_plot_dimensions_identical"] or not e["lighting_arms_match_exact_physical_lux_and_color"]:
        raise ValueError("Source same original physical site identity not verified")
    if e["independent_light_experiment_count"]!=1 or e["paper_publication_count"]!=2:
        raise ValueError("Treat two papers as one original lighting experiment")
    if e["source_published_plot_ID_codelist_crosspaper_verified"] or e["original_exact_plot_ID_and_effect_covariance_not_yet_verified"] is not True:
        raise ValueError("No original exact plot identity/SE covariance joining was performed")
    g=registry["new_joint_cross_trophic_estimate_gate"]
    if not g["shared_treatment_only_does_not_identify_spider_mediated_plant_phenology"]:
        raise ValueError("Source treatment arm cross-paper match is not causal mediation")
    if g["independent_crosspaper_same_plot_match_performed"]:
        raise ValueError("Original source plot/year data must be acquired before claiming a real plot join")
    if not g["herbivory_not_measured_by_Davies_pitfalls"]:
        raise ValueError("Source only sampled ground-active spiders and beetles")
    return {
        "schema":"alan_20261119_exeter_original_same_30_lightplots_two_trophic_outcomes_source_gate_v1",
        "field_programme":registry["programme_id"],
        "published_article_count":len(pubs),
        "independently_randomized_light_field_experiments":1,
        "independently_randomized_physical_field_plots":30,
        "physical_plots_per_lighting_regime":6,
        "light_arms":{
            x["canonical"]:{"Bennie2018":x["Bennie2018"],"Davies2017":x["Davies2017"],
                              "physical_plots":6}
            for x in arms
        },
        "physical_year_overlap":site["shared_field_years"],
        "original_plant_species_with_flowering_dates_n":3,
        "original_spider_beetle_caught_individuals":5180,
        "original_spider_beetle_taxa":136,
        "original_study_outcomes":{
            "Bennie2018":"grass flowering date, cover and biomass",
            "Davies2017":"day/night Araneae and Coleoptera pitfall assemblage catch"
        },
        "original_study_plants_and_spiders_same_light_plots":True,
        "original_plot_ID_codelist_usable_for_joint_source_model":False,
        "source_programme_can_be_used_as_two_meta_independent_studies":False,
        "actual_source_joint_cross_trophic_paired_model_fitted":False,
        "no_predator_mediation_or_species_seed_effect_estimated":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=DEFAULT)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--access",type=Path,required=False)
    a=p.parse_args()
    source=json.loads(a.source.read_text())
    result=validate(source)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("EXETER30_PROGRAMME_LIGHT_UNIT_LINK",result)
    if a.access:
        access=source_scan()
        a.access.parent.mkdir(parents=True,exist_ok=True)
        a.access.write_text(json.dumps(access,ensure_ascii=False,indent=2)+"\n")
        print("EXETER30_PUBLIC_ARCHIVE_SOURCE_ACCESS",access["public_archival_page_access"])


if __name__=="__main__":
    main()
