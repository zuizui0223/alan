#!/usr/bin/env python3
"""Macro audit across Page2021 original 623 pollinator taxa: visitation != function.

Only same SOURCE study, same plant, same source year/period, same SVE metric
and same visitation-rate units can form an individual plant functional
pollination network. This deliberately restricts SVE to author-reported
'number of conspecific pollen grains' and visit rates to stated comparable
per-unit-time, per-unit-area field efforts. q*V is a hypothetical effective
pollen delivery rate within a source group, not actual fruit/seed output.

Original Page et al. AJB DOI 10.1002/ajb2.1764
Figshare project 120921 article 16160265 file 30155367, sha locked.
No ALAN treatment, no nocturnality measurement, no plant-local extinction.
One published original study may contain multiple plant-years; cluster
source-study-level descriptive uncertainty rather than treating 623 taxa
as 623 independent experiments.

Nymphalidae Vanessa carye is incorrectly labelled 'moth' in 2 original
records. Do not overwrite original; run with and without this *verified*
butterfly mislabel. Harvard MCZ taxonomy primary source.
"""
from __future__ import annotations
import argparse, collections, json, math, random, statistics
from pathlib import Path
from scripts.analysis.alan_page2021_original_623_moth_fly_sameplant import (
    get_source,num,VKEY
)

SCHEMA="alan_page2021_original_pollinator_rank_vs_stigma_pollen_flux_macro_v1"
B=3999
SEED=20261118
VALID_VISITS={
 "number of visits per unit time per unit area",
 "number of flowers visited per unit time per unit area"
}
EFF="number of conspecific pollen grains"
UNITS="mean"

def candidate_groups(rows,taxon_validate=False):
    groups=collections.defaultdict(list)
    for r in rows:
        if (taxon_validate and r["pollinator.group"].strip()=="moth" and
            r["pollinator.taxon"].strip()=="Vanessa carye"):
            continue
        key=tuple(r[k].strip() for k in VKEY)
        groups[key].append(r)
    output=[]
    excluded=collections.Counter()
    for k,values in groups.items():
        if k[3].lower()!=EFF or k[4].lower()!=UNITS or k[5].lower() not in VALID_VISITS:
            excluded["not_conspecific_pollen_with_standard_visitation_rate"]+=1
            continue
        if len(values)<2:
            excluded["fewer_than_two_source_pollinator_rows"]+=1
            continue
        if len(set((r["pollinator.group"],r["pollinator.taxon"]) for r in values))!=len(values):
            excluded["duplicate_original_pollinator_taxa_same_source_cell"]+=1
            continue
        observed=[]
        for row in values:
            q=num(row["effectiveness.average"])
            rate=num(row["visit.value"])
            if q is None or rate is None or q<0 or rate<0:
                observed=[];break
            if rate>0:
                observed.append({
                  "taxon":row["pollinator.taxon"],
                  "guild":row["pollinator.group"],
                  "observed_visitation_rate":rate,
                  "source_single_visit_stigma_conspecific_pollen_mean":q,
                  "estimated_conspecific_pollen_delivery_rate":q*rate})
        if len(observed)<2:
            excluded["missing_numeric_or_at_least_two_active_visitors"]+=1
            continue
        totalV=sum(a["observed_visitation_rate"] for a in observed)
        totalP=sum(a["estimated_conspecific_pollen_delivery_rate"] for a in observed)
        if totalV<=0 or totalP<=0:
            excluded["zero_total_pollen_delivery_no_functional_rank"]+=1
            continue
        # Stable taxon name tie breaking is reported as a flagged ambiguity.
        orderedV=sorted(observed,key=lambda x:(-x["observed_visitation_rate"],x["taxon"]))
        orderedP=sorted(observed,key=lambda x:(-x["estimated_conspecific_pollen_delivery_rate"],x["taxon"]))
        tiesV=sum(abs(x["observed_visitation_rate"]-orderedV[0]["observed_visitation_rate"])<1e-12 for x in observed)
        tiesP=sum(abs(x["estimated_conspecific_pollen_delivery_rate"]-orderedP[0]["estimated_conspecific_pollen_delivery_rate"])<1e-12 for x in observed)
        if tiesV>1 or tiesP>1:
            excluded["ambiguous_equal_top_visitor_or_functional_taxon"]+=1
            continue
        rare=[{
            "original_taxon":a["taxon"],
            "original_guild":a["guild"],
            "visitation_share":a["observed_visitation_rate"]/totalV,
            "functional_pollen_rate_share":a["estimated_conspecific_pollen_delivery_rate"]/totalP
        } for a in observed if a["observed_visitation_rate"]/totalV<.05 and
                           a["estimated_conspecific_pollen_delivery_rate"]/totalP>=.10]
        topV=orderedV[0];topP=orderedP[0]
        def n_eff(z):
            w=[max(x[z],0) for x in observed]
            total=sum(w)
            return total**2/sum(a*a for a in w) if total>0 else None
        output.append({
            "source_study":k[0],"source_plant":k[1],"source_study_year_period":k[2],
            "original_SVE_measure":k[3],"original_visitation_rate_unit":k[5],
            "original_source_pollinator_taxa":len(observed),
            "observed_source_pollinator_guilds":sorted(set(a["guild"] for a in observed)),
            "top_visitor_taxon":topV["taxon"],"top_functional_taxon":topP["taxon"],
            "top_visitor_equals_top_functional_pollen":topV["taxon"]==topP["taxon"],
            "top_visitor_visit_rate_share":topV["observed_visitation_rate"]/totalV,
            "top_visitor_functional_pollen_share":
                 topV["estimated_conspecific_pollen_delivery_rate"]/totalP,
            "top_pollen_taxon_visit_share":topP["observed_visitation_rate"]/totalV,
            "top_pollen_taxon_function_share":
                 topP["estimated_conspecific_pollen_delivery_rate"]/totalP,
            "Hill2_effective_taxon_number_by_visit_rate":n_eff("observed_visitation_rate"),
            "Hill2_effective_taxon_number_by_pollen_delivery":n_eff("estimated_conspecific_pollen_delivery_rate"),
            "rare_under5pct_visitation_over10pct_functional_pollen":rare,
            "rare_functional_taxa_n":len(rare),
            "source_zero_single_visit_pollen_taxa_n":
                sum(a["source_single_visit_stigma_conspecific_pollen_mean"]==0
                    for a in observed),
            "top_visitor_source_pollen_efficiency_zero":
                topV["source_single_visit_stigma_conspecific_pollen_mean"]==0,
            "source_single_site_plants_not_independent_studies":True
        })
    return output,dict(excluded),len(groups)

def boot_study_rate(groups):
    ids=sorted(set(x["source_study"] for x in groups))
    if len(ids)<5:return {"status":"insufficient_source_independent_studies_STOP","source_studies":len(ids)}
    by=collections.defaultdict(list)
    for x in groups:by[x["source_study"]].append(x)
    rng=random.Random(SEED)
    rates=[]
    for _ in range(B):
        draw=[x for st in rng.choices(ids,k=len(ids)) for x in by[st]]
        rates.append(sum(not x["top_visitor_equals_top_functional_pollen"]
                         for x in draw)/len(draw))
    rates.sort()
    return {"status":"source_study_clustered_observational_descriptive_bootstrap",
       "source_unique_study_IDs":len(ids),"original_plant_groups":len(groups),
       "bootstrap_draws":B,
       "mismatch_fraction95":[rates[int(.025*(B-1))],rates[int(.975*(B-1))]],
       "does_not_estimate_ALAN_effect":True}

def run(rows):
    original,excluded,total_source_cells=candidate_groups(rows,False)
    taxon_validated,excluded_verified,total_verified_cells=candidate_groups(rows,True)
    def describe(x,exclusions,cells):
        m=[z for z in x if not z["top_visitor_equals_top_functional_pollen"]]
        rare=[z for z in x if z["rare_functional_taxa_n"]>0]
        functional_narrow=[z for z in x if z["Hill2_effective_taxon_number_by_pollen_delivery"]<
                                             z["Hill2_effective_taxon_number_by_visit_rate"]]
        bystudy=collections.defaultdict(list)
        for z in x:bystudy[z["source_study"]].append(z)
        studyrates=[sum(not a["top_visitor_equals_top_functional_pollen"] for a in v)/len(v)
                    for v in bystudy.values()]
        largest=sorted(bystudy.items(),key=lambda z:len(z[1]),reverse=True)
        other_only=[z for z in x if z["source_zero_single_visit_pollen_taxa_n"]==0]
        std_without=[sum(not a["top_visitor_equals_top_functional_pollen"] for a in x
                if a["source_study"]!=sid)/
                sum(a["source_study"]!=sid for a in x)
                for sid in bystudy if sum(a["source_study"]!=sid for a in x)>0]
        return {
          "original_source_study_plant_period_metric_visitation_cells":cells,
          "source_study_balanced_mean_mismatch_fraction":
               statistics.fmean(studyrates) if studyrates else None,
          "largest_contributing_source_study":{
              "original_study_id":largest[0][0],
              "plant_networks":len(largest[0][1]),
              "within_source_study_mismatch_fraction":
                  sum(not a["top_visitor_equals_top_functional_pollen"]
                      for a in largest[0][1])/len(largest[0][1])
          } if largest else None,
          "leave_any_single_original_study_out_mismatch_fraction_min_max":
              [min(std_without),max(std_without)] if std_without else None,
          "positive_q_only_visitor_taxa_groups_count":len(other_only),
          "positive_q_only_visitor_taxa_mismatch_n":sum(
              not a["top_visitor_equals_top_functional_pollen"]
              for a in other_only),
          "positive_q_only_visitor_taxa_mismatch_fraction":(
              sum(not a["top_visitor_equals_top_functional_pollen"]
                  for a in other_only)/len(other_only)) if other_only else None,
          "all_source_q_variances_or_visitation_rate_SEs_missing":True,
          "eligible_comparable_functional_pollen_plant_groups":len(x),
          "source_unique_studies":len({a["source_study"] for a in x}),
          "source_unique_plant_taxon_names":len({a["source_plant"] for a in x}),
          "top_visitor_differs_from_top_function_n":len(m),
          "top_visitor_differs_from_top_function_fraction":len(m)/len(x) if x else None,
          "rare_visitor_under5percent_but_over10percent_pollen_groups_n":len(rare),
          "source_pollen_network_narrower_than_visitation_network_n":
              len(functional_narrow),
          "original_groups_excluded_by_source_quality":exclusions,
          "source_group_study_unit_bootstrap":boot_study_rate(x),
          "examples_mismatched_functional_taxa":m[:12],
          "examples_rare_effectiveness_taxa":rare[:12],
          "all_eligible_source_plant_groups":x,
          "no_page_study_used_ALAN_lighting_intervention":True,
          "no_crossplant_pollen_effectiveness_transfer_allowed":True
        }
    return {
      "schema":SCHEMA,
      "original_Page_taxonlevel_source_n":len(rows),
      "original_Page_source_CSV_SHA256":"e28ebe99e0dfb27b73fbe6638e568bf98fa95b92c920e1d3d981e8eaea27519a",
      "original_source_group_only_pollen_average_and_standard_visit_units":True,
      "taxonomic_source_Vanessa_carye_Nymphalidae_butterfly_not_moth":True,
      "source_original_no_taxon_change":describe(original,excluded,total_source_cells),
      "excluding_independently_verified_Vanessa_mislabeled_moth":describe(
          taxon_validated,excluded_verified,total_verified_cells),
      "source_not_all_pollinator_taxa_independently_verified":True,
      "source_no_pollen_SDs_or_ALAN_treatments":True,
      "source_original_study_year_label_can_cover_multiple_calendar_years":True,
      "source_most_important_caveat":"q times V is an observed SVE x visit-rate proxy, not seed set or actual ALAN response"
    }
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    q=p.parse_args()
    data=run(get_source())
    f=Path(q.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    for title,key in (("raw","source_original_no_taxon_change"),
                      ("verified","excluding_independently_verified_Vanessa_mislabeled_moth")):
        x=data[key]
        print("PAGE2021_FUNCTIONAL_VISITOR_MISMATCH_SOURCE",json.dumps({
          "version":title,"source_plant_groups":x["eligible_comparable_functional_pollen_plant_groups"],
          "independent_studies":x["source_unique_studies"],
          "plants":x["source_unique_plant_taxon_names"],
          "mismatch_n":x["top_visitor_differs_from_top_function_n"],
          "mismatch_fraction":x["top_visitor_differs_from_top_function_fraction"],
          "rare_important_groups":x["rare_visitor_under5percent_but_over10percent_pollen_groups_n"],
          "study_balanced_mismatch":x["source_study_balanced_mean_mismatch_fraction"],
          "largest_study":x["largest_contributing_source_study"],
          "leave_one_study_range":x["leave_any_single_original_study_out_mismatch_fraction_min_max"],
          "positive_q_only_n":x["positive_q_only_visitor_taxa_groups_count"],
          "positive_q_only_mismatch":x["positive_q_only_visitor_taxa_mismatch_n"],
          "positive_q_only_rate":x["positive_q_only_visitor_taxa_mismatch_fraction"],
          "source_exclusions":x["original_groups_excluded_by_source_quality"],
          "study_bootstrap":x["source_group_study_unit_bootstrap"],
          "top_mismatch_examples":x["examples_mismatched_functional_taxa"][:8],
          "rare_important_examples":x["examples_rare_effectiveness_taxa"][:8]
        },ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
