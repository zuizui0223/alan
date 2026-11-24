#!/usr/bin/env python3
"""Narrow cross-programme taxonomic-composition gate, not plant mediation.

Pair two physical experimental programmes that both measure aquatic
invertebrate drift, but restrict *numerator* to Ephemeroptera+Trichoptera and
*denominator* to Diptera+Coleoptera. Other insects and non-insect animals
are explicitly excluded, never called zeros.

Henn 2012: 5 natural rivers, one dark and one lit sampling night.
EERES 2021: 6 shared header tanks, repeated flume samples.
The ratio of original arm means across independent blocks is the common
estimand; mean within-block log odds is a separate Henn sensitivity.
This is a source-gated descriptive comparison, NOT a k=2 accepted meta.
"""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import statistics
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[2]
HENN=ROOT/"database/alan_experiments/arthropod_expansion_20261122/stream_drift_20261124"
EERES_DOI="10.5061/dryad.rfj6q57k5"
XLSX_SHA="496dcfa07436619084dc677317525dc674775e9965a0fdfb482cfa1cb4d8977f"
NAME="DISCONNECTED_alldata.xlsx"
source_path=ROOT/"scripts/acquisition/alan_eeres_2021_source_gate.py"
spec=importlib.util.spec_from_file_location("eeres_source_gate",source_path)
source_mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(source_mod)

# Curated source-header-to-order classifications. All source headers must match.
MAYFLY={"Baetis","Ecdyonurus","Ephemerella","Leptophlebiidae"}
CADDIS={"Brachicentridae","Hydropsyche","Limnephilidae","Odontocerum",
        "Philopotamus","Rhyacophila","Sericostoma"}
DIPTERA={"Ceratopogonidae","Chironomini","Dixidae","Empididae","Limoniidae",
         "Orthocladiinae","Psychodidae","Simuliidae","Stratiomyidae",
         "Tanypodinae","Tanytarsini"}
COLEOPTERA={"Elmis","Elodes","Hydraena","Hydroporus","Laccophilus","Limnius","Oreodytes"}
STONEFLY={"Leuctra","Perlodes","Protonemura"}
NON_INSECT={"Gammarus","Hydrachnidia","Lumbricidae","Mermithidae","Naididae","Polycelis"}
TAXON_GROUPS={
    "E_T":MAYFLY|CADDIS,
    "D_C":DIPTERA|COLEOPTERA,
    "Plecoptera_excluded":STONEFLY,
    "NonInsect_excluded":NON_INSECT,
}
TCRIT={5:2.7764451051977987,6:2.570581835636305}

def source_taxon_mapping(source_headers):
    if len(source_headers)!=38 or len(set(source_headers))!=38:
        raise ValueError("STOP: expected 38 unique source original drift headers")
    known=set.union(*TAXON_GROUPS.values())
    if source_headers!=set(source_headers):
        pass  # source_headers may be a sequence, validate sets below
    if known!=set(source_headers):
        raise ValueError(f"STOP: source taxonomy changed: missing={known-set(source_headers)}, new={set(source_headers)-known}")
    if sum(len(v) for v in TAXON_GROUPS.values())!=len(known):
        raise ValueError("STOP: ambiguous multi-assigned taxon")
    return {k:g for g,names in TAXON_GROUPS.items() for k in names}

def ratio_from_block_arms(block_rows,programme):
    """A vector of per-independent-block (A_lit,B_lit,A_dark,B_dark)."""
    if len(block_rows) not in TCRIT:
        raise ValueError("STOP: expected n=5 river or n=6 header independent units")
    n=len(block_rows)
    for a,b,c,d in block_rows:
        if any(not math.isfinite(v) or v<0 for v in (a,b,c,d)):
            raise ValueError("STOP: invalid negative, nonfinite or missing block counts")
    mu=[statistics.mean(row[j] for row in block_rows) for j in range(4)]
    if min(mu)<=0:
        raise ValueError("STOP: group mean zero, cannot form log ratio")
    logrr=math.log(mu[0])-math.log(mu[1])-math.log(mu[2])+math.log(mu[3])
    grad=[1/mu[0],-1/mu[1],-1/mu[2],1/mu[3]]
    scores=[sum(g*(r[j]-mu[j]) for j,g in enumerate(grad)) for r in block_rows]
    se=statistics.stdev(scores)/math.sqrt(n)
    critical=TCRIT[n]
    return {
        "programme":programme,"n_independent_physical_units":n,
        "arm_means_A_lit_B_lit_A_dark_B_dark":mu,
        "four_order_log_composition_RRR":logrr,
        "four_order_composition_RRR":math.exp(logrr),
        "delta_SE_shared_physical_units":se,
        "descriptive_t95":[math.exp(logrr-critical*se),math.exp(logrr+critical*se)],
        "block_unit_data_is_original_source_sample_averages_not_insect_individual_replication":True,
        "block_marginals":[{"A_lit":a,"B_lit":b,"A_dark":c,"B_dark":d} for a,b,c,d in block_rows],
    }

def henn_study():
    with (HENN/"henn2014_derived_order_by_stream.csv").open(encoding="utf-8",newline="") as f:
        rows={x["order"]:x for x in csv.DictReader(f)}
    if set(rows)!={"Coleoptera","Diptera","Ephemeroptera","Lepidoptera","Odonata","Trichoptera"}:
        raise ValueError("STOP: original Henn source order list changed")
    river_ids=("Barton","Bull","Onion","Comal","SanMarcos")
    expected_lit=(57,49,46,286,793)
    expected_dark=(7,59,98,601,1196)
    matched=[]
    for i,r in enumerate(river_ids):
        a=sum(int(rows[o][r+"_lit"]) for o in ("Ephemeroptera","Trichoptera"))
        b=sum(int(rows[o][r+"_lit"]) for o in ("Coleoptera","Diptera"))
        c=sum(int(rows[o][r+"_dark"]) for o in ("Ephemeroptera","Trichoptera"))
        d=sum(int(rows[o][r+"_dark"]) for o in ("Coleoptera","Diptera"))
        total_lit=sum(int(rows[o][r+"_lit"]) for o in rows)
        total_dark=sum(int(rows[o][r+"_dark"]) for o in rows)
        if (total_lit,total_dark)!=(expected_lit[i],expected_dark[i]):
            raise ValueError("STOP: Henn family-derived order sums changed")
        matched.append((a,b,c,d))
    result=ratio_from_block_arms(matched,"Henn2014_TexasFive_natural_streams")
    if any(min(row)<=0 for row in matched):
        raise ValueError("Henn source 4-order counts should all be positive")
    block_logs=[math.log(a/b)-math.log(c/d) for a,b,c,d in matched]
    mu=statistics.mean(block_logs)
    se=statistics.stdev(block_logs)/math.sqrt(5)
    result.update({
        "alternate_equal_stream_geometric_RRR":math.exp(mu),
        "alternate_equal_stream_descriptive_t95":[math.exp(mu-TCRIT[5]*se),math.exp(mu+TCRIT[5]*se)],
        "streams_with_four_order_RRR_above_one":sum(r>0 for r in block_logs),
        "table_total_discrepancy_lit_count":2,
        "source_type":"published_original_Table_1_manual_derived_orders",
    })
    return result

def source_original_xlsx():
    """Use the existing source gate's official API+redirect allowlist and original hash."""
    src=source_mod.run()
    if src.get("status")!="SOURCE_WORKBOOK_VERIFIED":
        raise ValueError("STOP: genuine EERES Dryad workbook not accessible")
    if src.get("original_workbook",{}).get("sha256")!=XLSX_SHA:
        raise ValueError("STOP: source workbook hash changed")
    versions=[u for u in src.get("metadata_links",[])
              if u.startswith("https://datadryad.org/api/v2/versions/") and u.rstrip("/").split("/")[-1].isdigit()]
    if not versions:
        raise ValueError("STOP: source dataset metadata has no official version link")
    for u in versions:
        try:
            original_zip=source_mod.grab(u+"/download")
            with zipfile.ZipFile(io.BytesIO(original_zip)) as arch:
                candidates=[x for x in arch.namelist() if x.rsplit("/",1)[-1]==NAME]
                if len(candidates)!=1:
                    continue
                raw=arch.read(candidates[0])
                if hashlib.sha256(raw).hexdigest()==XLSX_SHA:
                    return raw
        except (OSError,ValueError,zipfile.BadZipFile):
            continue
    raise ValueError("STOP: original signed XLSX missing from official source")

def eeres_study(original):
    import openpyxl
    wb=openpyxl.load_workbook(io.BytesIO(original),read_only=True,data_only=True)
    ws=wb["drift"];rows=ws.iter_rows(values_only=True)
    header=[str(x or "").strip() for x in next(rows)]
    if header[1:3]!=["HT","flume"] or "treatment_up" not in header:
        raise ValueError("STOP: drift source structural columns changed")
    if len(header[:8])!=8:
        raise ValueError("STOP: wrong source taxon offset")
    names=header[8:]
    mapping=source_taxon_mapping(names)
    idx={k:header.index(k) for k in ("HT","flume","treatment_up")}
    column_groups={g:[8+j for j,n in enumerate(names) if mapping[n]==g] for g in TAXON_GROUPS}
    expected={"CT":12,"L":4,"F":4,"L+F":4}
    values={}
    flumes=set();nr=0
    for row in rows:
        if row[0] is None:continue
        ht=str(row[idx["HT"]]).strip()
        flume=str(row[idx["flume"]]).strip()
        tr=str(row[idx["treatment_up"]]).strip()
        if tr=="CT_up":tr="CT"
        if ht not in {f"ht{i}" for i in range(1,7)} or tr not in expected:
            raise ValueError("STOP: block/treatment code unexpected")
        if len(row)!=len(header):
            raise ValueError("STOP: taxon source record width changed")
        numbers=row[8:]
        if any(v is None or isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0 for v in numbers):
            raise ValueError("STOP: nonnumeric or missing actual source cells")
        sums={g:sum(row[j] for j in indices) for g,indices in column_groups.items()}
        values.setdefault((ht,tr),[]).append(sums)
        flumes.add(flume);nr+=1
    if nr!=144 or len(flumes)!=18:
        raise ValueError(f"STOP: EERES source grain changed rows={nr},flumes={len(flumes)}")
    blocks=[]
    exclusions={"plecoptera_counts_by_arm":{},"noninsect_counts_by_arm":{}}
    for i in range(1,7):
        ht=f"ht{i}"
        for tr,n in expected.items():
            if len(values.get((ht,tr),[]))!=n:
                raise ValueError(f"STOP: drift unbalanced source {ht}/{tr} expected={n}")
        def avg(tr,group):
            return statistics.mean(x[group] for x in values[ht,tr])
        blocks.append((avg("L","E_T"),avg("L","D_C"),avg("CT","E_T"),avg("CT","D_C")))
    result=ratio_from_block_arms(blocks,"EERES2021_six_header_tanks_controlled_flumes")
    # Fully factorial same-HT sensitivity: compare light effect on ET/DC drift
    # at normal flow (L/CT) with its effect under reduced flow (L+F/F).
    # All four treatments share each of the exact six physical header tanks.
    flow_blocks=[]
    for i in range(1,7):
        ht=f"ht{i}"
        def flow_avg(tr,group):
            return statistics.mean(x[group] for x in values[ht,tr])
        flow_blocks.append((
            flow_avg("L+F","E_T"), flow_avg("L+F","D_C"),
            flow_avg("F","E_T"), flow_avg("F","D_C"),
        ))
    flow_effect=ratio_from_block_arms(flow_blocks,
                "EERES2021_six_header_tanks_reduced_flow")
    base_arr=[(r["A_lit"],r["B_lit"],r["A_dark"],r["B_dark"]) for r in result["block_marginals"]]
    normal_mu=result["arm_means_A_lit_B_lit_A_dark_B_dark"]
    flow_mu=flow_effect["arm_means_A_lit_B_lit_A_dark_B_dark"]
    def influence(r,mu):
        g=(1/mu[0],-1/mu[1],-1/mu[2],1/mu[3])
        return sum(gg*(rr-mm) for gg,rr,mm in zip(g,r,mu))
    d_scores=[influence(flow_blocks[i],flow_mu)-influence(base_arr[i],normal_mu)
              for i in range(6)]
    d_se=statistics.stdev(d_scores)/math.sqrt(6)
    d=flow_effect["four_order_log_composition_RRR"]-result["four_order_log_composition_RRR"]
    result["flow_moderator_exploratory"]={
        "contrast":"log[ (ET/DC)_LF / (ET/DC)_F ] - log[ (ET/DC)_L / (ET/DC)_CT ]",
        "reduced_flow_light_composition_RRR":flow_effect["four_order_composition_RRR"],
        "reduced_flow_logRRR":flow_effect["four_order_log_composition_RRR"],
        "reduced_flow_delta_SE":flow_effect["delta_SE_shared_physical_units"],
        "reduced_flow_descriptive_t5_95":flow_effect["descriptive_t95"],
        "within_same_six_HT_moderation_logRRR":d,
        "moderation_ratio_of_RRRs":math.exp(d),
        "joint_shared_HT_delta_SE":d_se,
        "descriptive_t5_95_moderator":[math.exp(d-2.570581835636305*d_se),
                                      math.exp(d+2.570581835636305*d_se)],
        "n_true_independent_header_tanks":6,
        "not_a_water_wet_vs_dry_gradient_or_plant_fitness":True,
        "posthoc_not_confirmatory":True,
    }
    result.update({
        "source_sha256":XLSX_SHA,
        "source_drift_taxon_columns":len(names),
        "source_excluded_groups_counts":{"Plecoptera":len(STONEFLY),"non_insect":len(NON_INSECT)},
        "included_group_taxon_columns":{"E_T":len(TAXON_GROUPS["E_T"]),"D_C":len(TAXON_GROUPS["D_C"])},
        "physical_header_n":6,"original_source_rows":nr,
        "source_type":"Dryad_original_xlsx_SHA_verified",
    })
    return result

def analyze(eeres_original):
    a=henn_study()
    b=eeres_study(eeres_original)
    return {
        "comparison":"Two separate true physical programmes, in-stream aquatic insect drift composition E+T relative to Diptera+Coleoptera",
        "status":"SOURCE_AUTHENTICATED_EXPLORATORY_TWO_PROGRAMME_COMPARISON_NOT_ACCEPTED_META",
        "meta_k_added":0,"mkh_or_reml_pooled_not_fitted":True,
        "source_dois":["10.1080/02705060.2014.900654",EERES_DOI],
        "henn":a,"eeres":b,
        "warning":"Both are post hoc; Henn uses different nights and extreme incandescent illumination in five natural rivers, EERES uses six shared experimental flume headers, multiple weeks, moderate ALAN. Source Henn family sum differs from published Total by +2 for San Marcos lit. Taxon labels differ in resolution; excluded three Plecoptera and six noninsects in EERES for shared-four-order comparison. Comparing signs across two experiments is not a meta causal mechanism nor a wet-vs-dry or plant-fitness moderator."
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=Path("results/alan_henn_eeres_four_order_drift_20261124.json"))
    args=ap.parse_args()
    out=analyze(source_original_xlsx())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("SOURCE_MATCHED_FOUR_ORDER",json.dumps({
        "henn":{k:out["henn"][k] for k in ("four_order_composition_RRR","four_order_log_composition_RRR","delta_SE_shared_physical_units","descriptive_t95","n_independent_physical_units","alternate_equal_stream_geometric_RRR")},
        "eeres":{k:out["eeres"][k] for k in ("four_order_composition_RRR","four_order_log_composition_RRR","delta_SE_shared_physical_units","descriptive_t95","n_independent_physical_units","flow_moderator_exploratory")},
        "meta_k_added":out["meta_k_added"],
    },ensure_ascii=False))

if __name__=="__main__":
    main()
