#!/usr/bin/env python3
"""Family-level explanation of *other* catches in Sullivan aquatic emergence traps.

Source-author original 2013 wetland file, SHA-pinned by the existing gate.
This is an exploratory description of SOURCE CATCH COMPOSITION by original
habitat block, not proof that every non-Chironomidae adult emerged from water.
The source AquaticEmergent sample also caught taxa with terrestrial immature
stages. A contrasting response of the broad complement could be due to
different capture/attraction or habitat, not a wetland predator mechanism.

No raw sample-level author file is redistributed.
"""
from __future__ import annotations
import argparse,csv,io,importlib.util,json,math,itertools,statistics
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
M=ROOT/"scripts/analysis/alan_sullivan2019_chironomid_source_gate_20261124.py"
spec=importlib.util.spec_from_file_location("sullivan_original_family",M)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
LIGHTS=("LEDI","LEDH","HPSI","HPSH")

def inspect(original):
    if base.DIGEST!="385444ab318ec4b15edd3161c740aa7f4b2d9b973fe4c2ddf6c29cc56cdec484":
        raise ValueError("STOP source original SHA changed")
    src=list(csv.DictReader(io.StringIO(original.decode("utf-8-sig")),delimiter="\t"))
    if len(src)!=1161:raise ValueError("STOP original author rows not 1161")
    aq=[r for r in src if r["Sample"].strip()=="AquaticEmergent"]
    if len(aq)!=198:raise ValueError("STOP original aquatic records 198")
    source_counts=Counter()
    blocks=defaultdict(Counter)
    for r in aq:
        group,tr,fam=r["Group"].strip(),r["Treatment"].strip(),r["Taxa"].strip()
        if group not in ("1","2","3") or tr not in ("DRK",)+LIGHTS or not fam:
            raise ValueError("STOP original block/treatment/family invalid")
        x=float(r["Count"])
        if not math.isfinite(x) or x<0 or not x.is_integer():
            raise ValueError("STOP source taxon-count nonnegative integer required")
        n=int(x)
        source_counts[fam]+=n
        blocks[group,tr][fam]+=n
    if source_counts["Chironomidae"]!=4960 or sum(source_counts.values())!=5372:
        raise ValueError("STOP original source family total drift")
    for g in ("1","2","3"):
        for tr in ("DRK",)+LIGHTS:
            if (g,tr) not in blocks:
                raise ValueError("STOP source physical wetland sampling cell missing")

    nonchiro=[(fam,n) for fam,n in source_counts.most_common() if fam!="Chironomidae"]
    main_other=nonchiro[:6]
    chosen=[fam for fam,_ in main_other]
    per_block=[]
    chiro_vs_mosquito_logs=[]
    chiro_vs_mosquito_blocks=[]
    for group in ("1","2","3"):
        d=blocks[group,"DRK"]
        l=Counter()
        for tr in LIGHTS:l.update(blocks[group,tr])
        d_total=sum(d.values());l_total=sum(l.values())
        d_ch=d["Chironomidae"];l_ch=l["Chironomidae"]
        if not (d_total>d_ch>0 and l_total>l_ch>0):
            raise ValueError("STOP physical block missing complementary taxa")
        # Prespecified only after inspecting the complement: Culicidae and
        # Chironomidae both have aquatic larval life histories, unlike the
        # heterogeneous catch remainder. This is POST HOC and has n=3 only.
        if not all(v>0 for v in (
            d["Chironomidae"],l["Chironomidae"],d["Culicidae"],l["Culicidae"]
        )):
            raise ValueError("STOP: zero catch within Chironomidae/Culicidae source family pair")
        cm_rrr=(l["Chironomidae"]/d["Chironomidae"])/(l["Culicidae"]/d["Culicidae"])
        chiro_vs_mosquito_logs.append(math.log(cm_rrr))
        chiro_vs_mosquito_blocks.append({
            "block":group,
            "source_chironomidae_dark":d["Chironomidae"],
            "source_chironomidae_four_lamps":l["Chironomidae"],
            "source_culicidae_dark":d["Culicidae"],
            "source_culicidae_four_lamps":l["Culicidae"],
            "two_family_RRR":cm_rrr,
        })
        per_block.append({
            "wetland_block":group,
            "dark_all_taxa_captured":d_total,
            "four_illumination_lamp_all_taxa_captured":l_total,
            "dark_chironomidae_captured":d_ch,
            "four_illumination_lamp_chironomidae_captured":l_ch,
            "dark_all_other_captured":d_total-d_ch,
            "four_illumination_lamp_all_other_captured":l_total-l_ch,
            "major_other_labels":[
                {"source_family_label":fam,"dark":d[fam],"four_lamp_total":l[fam],
                 "delta_proportion_of_all_other_light_minus_dark":(
                    l[fam]/(l_total-l_ch)-d[fam]/(d_total-d_ch))}
                for fam in chosen
            ],
            "all_other_total_by_lamp":{
                tr:sum(blocks[group,tr].values())-blocks[group,tr]["Chironomidae"] for tr in LIGHTS
            },
            "chiro_vs_all_other_response_RRR":(l_ch/d_ch)/((l_total-l_ch)/(d_total-d_ch)),
        })
    m=statistics.mean(chiro_vs_mosquito_logs)
    se=statistics.stdev(chiro_vs_mosquito_logs)/math.sqrt(3)
    t2=4.302652729911275
    signs=sum(1 for choice in itertools.product((-1,1),repeat=3)
              if abs(sum(x*y for x,y in zip(choice,chiro_vs_mosquito_logs))/3)>=abs(m)-1e-12)
    pair={
        "source_selected_families":["Chironomidae","Culicidae"],
        "source_selection_posthoc":True,
        "both_taxa_have_aquatic_immature_development":True,
        "three_original_block_counts_and_ratios":chiro_vs_mosquito_blocks,
        "physical_habitat_block_n":3,
        "equal_block_logRRR":m,
        "equal_block_logRRR_SE":se,
        "equal_block_RRR":math.exp(m),
        "descriptive_t2_95_RRR":[math.exp(m-t2*se),math.exp(m+t2*se)],
        "three_block_signflip_reference_p":signs/8,
        "not_independent_experimental_meta_effect":True,
        "source_lamp_menu_four_depend_on_same_dark_each_block":True,
        "multiple_candidate_family_comparisons_not_predeclared":True,
    }
    return {
        "source_original_TSV_sha256":base.DIGEST,
        "publication_doi":"10.1002/eap.1821",
        "source_data_doi":"10.5061/dryad.gb9534p",
        "true_wetland_habitat_blocks":3,
        "original_aquatic_emergence_catches":5372,
        "source_family_label_richness":len(source_counts),
        "source_total_non_chironomidae_catches":412,
        "source_non_chironomidae_family_counts":dict(nonchiro),
        "top_six_other_families":dict(main_other),
        "physical_block_family_contributions":per_block,
        "explicit_two_aquatic_larval_family_sensitivity":pair,
        "taxonomic_classification_limit":(
            "Non-Chironomidae source catch includes families with terrestrial immature development, "
            "not exclusively freshwater emerging insects. Family identity from author label, not individual voucher re-identification."
        ),
        "lighting_stations_equal_sampling_effort_confirmation_pending":True,
        "source_spatial_count_unit":"3 physical habitat blocks, not 15 independently illuminated habitat programmes",
        "new_meta_programmes_added":0,
        "do_not_apply_posthoc_family_selection_as_new_confirmation":True,
        "status":"ORIGINAL_FAMILY_DIAGNOSTIC_DESCRIPTIVE_ONLY",
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path,default=Path("results/alan_sullivan2019_other_emergence_family_drivers_20261124.json"))
    args=a.parse_args()
    x=inspect(base.get_original())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(x,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("SULLIVAN_OTHER_FAMILY_SOURCE",json.dumps({
      "data_authenticated":True,"other_taxa":x["source_non_chironomidae_family_counts"],
      "true_block_decomposition":x["physical_block_family_contributions"],
      "aquatic_larval_two_family_comparison":x["explicit_two_aquatic_larval_family_sensitivity"],
      "meta_k_increment":0
    },ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
