#!/usr/bin/env python3
"""Order-deletion fragility diagnostic for the Henn-EERES 2-programme drift meta.

Excludes each of the four originally matched orders in turn from E+T vs D+C.
Every exclusion CHANGES the estimand and is a post hoc diagnostic, not 4
new, independent hypotheses or a new effect-size study. Both original studies
are kept as one physical programme each, with the original 5 rivers and 6
header tanks as respective inference units. No taxon record redistribution.
"""
import argparse
import csv
import io
import json
import math
import statistics
import sys
from pathlib import Path
import importlib.util

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/"database/alan_experiments/arthropod_expansion_20261122/stream_drift_20261124"
SOURCE=ROOT/"scripts/analysis/alan_20261124_henn_eeres_taxon_matched_drift.py"
spec=importlib.util.spec_from_file_location("henn_eeres_base", SOURCE)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
ORDERS=("Ephemeroptera","Trichoptera","Diptera","Coleoptera")
ET=ORDERS[:2];DC=ORDERS[2:]
T1=12.706204736432095

def Henn():
    with (HERE/"henn2014_derived_order_by_stream.csv").open(encoding="utf-8",newline="") as f:
        items=list(csv.DictReader(f))
    if len(items)!=6 or set(x["order"] for x in items)!={"Coleoptera","Diptera","Ephemeroptera","Lepidoptera","Odonata","Trichoptera"}:
        raise ValueError("STOP: Henn source-order rows changed")
    table={x["order"]:x for x in items}
    original=base.henn_study()  # asserts all five six-order total rows
    ids=("Barton","Bull","Onion","Comal","SanMarcos")
    return [
        {o:(int(table[o][r+"_lit"]),int(table[o][r+"_dark"])) for o in ORDERS}
        for r in ids
    ],original

def EERES(author_xlsx):
    import openpyxl
    wb=openpyxl.load_workbook(io.BytesIO(author_xlsx),read_only=True,data_only=True)
    rows=wb["drift"].iter_rows(values_only=True)
    head=[str(x or "").strip() for x in next(rows)]
    if head[1:3]!=["HT","flume"] or "treatment_up" not in head or len(head)!=46:
        raise ValueError("STOP: original EERES drift schema mismatch")
    names=head[8:]
    base.source_taxon_mapping(names)
    taxo={"Ephemeroptera":base.MAYFLY,"Trichoptera":base.CADDIS,
          "Diptera":base.DIPTERA,"Coleoptera":base.COLEOPTERA}
    idx={o:[8+j for j,name in enumerate(names) if name in label] for o,label in taxo.items()}
    if {o:len(z) for o,z in idx.items()}!={"Ephemeroptera":4,"Trichoptera":7,"Diptera":11,"Coleoptera":7}:
        raise ValueError("STOP: changed original taxon-column counts")
    hti=head.index("HT");tri=head.index("treatment_up");fi=head.index("flume")
    vals={};flumes=set();nr=0
    for row in rows:
        if row[0] is None: continue
        ht=str(row[hti]).strip();tr=str(row[tri]).strip()
        if tr=="CT_up":tr="CT"
        if ht not in {f"ht{i}" for i in range(1,7)} or tr not in {"L","CT","F","L+F"}:
            raise ValueError("STOP: experiment label changed")
        if len(row)!=len(head) or any(not isinstance(v,(int,float)) or isinstance(v,bool) or not math.isfinite(v) or v<0 for v in row[8:]):
            raise ValueError("STOP: numeric source change")
        values={o:sum(row[j] for j in ind) for o,ind in idx.items()}
        vals.setdefault((ht,tr),[]).append(values)
        flumes.add(str(row[fi]).strip());nr+=1
    if nr!=144 or len(flumes)!=18:
        raise ValueError("STOP: physical sample count changed")
    result=[]
    for i in range(1,7):
        ht="ht"+str(i);pair={}
        for tr,expect in (("L",4),("CT",12),("F",4),("L+F",4)):
            if len(vals.get((ht,tr),[]))!=expect:
                raise ValueError(f"STOP: physical replicated sample count changed {ht} {tr}")
        for o in ORDERS:
            pair[o]=(statistics.mean(x[o] for x in vals[ht,"L"]),
                     statistics.mean(x[o] for x in vals[ht,"CT"]))
        result.append(pair)
    original=base.eeres_study(author_xlsx)
    return result,original

def build_arms(blocks,exclude):
    if exclude not in (None,*ORDERS):
        raise ValueError("STOP: unknown order deletion")
    num=set(ET)-({exclude} if exclude else set())
    den=set(DC)-({exclude} if exclude else set())
    if not num or not den:
        raise ValueError("STOP: empty numerator or denominator")
    arms=[]
    for r in blocks:
        arms.append((
            sum(r[o][0] for o in num),sum(r[o][0] for o in den),
            sum(r[o][1] for o in num),sum(r[o][1] for o in den)))
    return arms

def mKH_two(a,b):
    yi=[a["four_order_log_composition_RRR"],b["four_order_log_composition_RRR"]]
    vi=[a["delta_SE_shared_physical_units"]**2,b["delta_SE_shared_physical_units"]**2]
    diff=yi[0]-yi[1]
    tau2=max(0.,(diff*diff-vi[0]-vi[1])/2)
    w=[1/(z+tau2) for z in vi]
    mu=sum(yy*ww for yy,ww in zip(yi,w))/sum(w)
    q=sum(ww*(yy-mu)**2 for yy,ww in zip(yi,w))
    se=math.sqrt(max(1.,q)/sum(w))
    return {"k_independent_programmes":2,
            "pooled_RRR":math.exp(mu),
            "mKH_t1_95":[math.exp(mu-T1*se),math.exp(mu+T1*se)],
            "reml_tau2":tau2,
            "posthoc_not_multiple_independent_confirmations":True}


def independent_block_order_ratios(blocks):
    """Absolute capture/counted drift responses for interpretation only.

    These are *not* flow-standardized between studies, not functional fluxes,
    and not evidence of changing true populations. n is genuine block count.
    """
    n=len(blocks)
    critical={5:2.7764451051977987,6:2.570581835636305}.get(n)
    if critical is None:raise ValueError("STOP: unknown physical block count")
    out={}
    for order in ORDERS:
        light=[float(row[order][0]) for row in blocks]
        dark=[float(row[order][1]) for row in blocks]
        ml=statistics.mean(light);md=statistics.mean(dark)
        if ml<=0 or md<=0:raise ValueError("STOP: cannot infer order RR with no source counts")
        logratio=math.log(ml/md)
        influences=[(a-ml)/ml-(b-md)/md for a,b in zip(light,dark)]
        se=statistics.stdev(influences)/math.sqrt(n)
        out[order]={"lit_over_dark_catch_RR":math.exp(logratio),
                    "source_mean_lit":ml,"source_mean_dark":md,
                    "n_independent_physical_blocks":n,
                    "descriptive_t95":[math.exp(logratio-critical*se),
                                       math.exp(logratio+critical*se)],
                    "paired_block_delta_se":se}
    return out

def run(original_xlsx):
    hh,hbase=Henn()
    ee,ebase=EERES(original_xlsx)
    values=[]
    for omitted in (None,*ORDERS):
        hs=base.ratio_from_block_arms(build_arms(hh,omitted),"Henn2012_5_streams")
        es=base.ratio_from_block_arms(build_arms(ee,omitted),"EERES2021_6_header_tanks")
        values.append({
            "deleted_order":omitted or "none",
            "comparison":"+".join(o for o in ET if o!=omitted)+"/"+"+".join(o for o in DC if o!=omitted),
            "henn":{"RRR":hs["four_order_composition_RRR"],"paired_se":hs["delta_SE_shared_physical_units"],
                    "descriptive_t95":hs["descriptive_t95"],"physical_n":5},
            "eeres":{"RRR":es["four_order_composition_RRR"],"paired_se":es["delta_SE_shared_physical_units"],
                     "descriptive_t95":es["descriptive_t95"],"physical_n":6},
            "exploratory_two_programme_mKH":mKH_two(hs,es),
        })
        if omitted is None:
            for observed,ref in ((hs,hbase),(es,ebase)):
                if not math.isclose(observed["four_order_composition_RRR"],ref["four_order_composition_RRR"],rel_tol=1e-12):
                    raise ValueError("STOP: primary main four-order RRR was not reproduced")
                if not math.isclose(observed["delta_SE_shared_physical_units"],ref["delta_SE_shared_physical_units"],rel_tol=1e-12):
                    raise ValueError("STOP: original primary sampling variance not reproduced")
            if not math.isclose(values[-1]["exploratory_two_programme_mKH"]["pooled_RRR"],1.2159960467857396,rel_tol=1e-12):
                raise ValueError("STOP: original 1.216 meta result not reproduced")
    return {
        "analysis":"leave-one-insect-order-out posterior posthoc sensitivity, NOT a new systematic study",
        "source_henn":"Original Henn 2014 Table 1 hand-derived orders,5 paired rivers;2-individual source discrepancy retained",
        "source_eeres_sha256":base.XLSX_SHA,
        "physical_k":[5,6],"new_independent_studies":0,
        "per_order_catch_ratios_not_drift_density":{"henn":independent_block_order_ratios(hh),
                "eeres":independent_block_order_ratios(ee)},
        "original_meta_unchanged":True,
        "source_taxonomic_grain":"Henn source is family aggregate by order; EERES source genus/family/subfamily/tribe labels classified into order",
        "comparisons":values,
        "limitations":[
          "Four posthoc deletions change the relative composition estimand; no multiple-testing corrected confirmation",
          "A low base count may make t and delta approximate SE unreliable for some deletions",
          "EERES source taxonomy is incomplete by species; Henn family Table 1 source has original totals discrepancies",
          "Henn high-illuminance one-night natural stream vs EERES lower-illuminance weeks-long experimental flume",
          "These results are insect drift composition, not overall drift abundance, cross-boundary flux, fish predation or wetland plant reproductive fitness"
        ],
    }


def verify_frozen(result, frozen):
    """Fail closed on any drift in the real author-source derived five comparisons."""
    if frozen.get("eeres_source_sha256")!=result["source_eeres_sha256"]:
        raise ValueError("STOP: original Dryad XLSX fingerprint mismatch")
    if frozen.get("independent_programmes")!=2 or not frozen.get("posthoc_deletions_change_estimand"):
        raise ValueError("STOP: invalid effect identity or independent study count")
    expected=frozen.get("rows",[])
    observed=result.get("comparisons",[])
    if len(expected)!=5 or len(observed)!=5:
        raise ValueError("STOP: missing one original order-deletion comparison")
    for a,b in zip(expected,observed):
        if a["deleted_order"]!=b["deleted_order"]:
            raise ValueError("STOP: original order-deletion identity drift")
        for key,val in (
            ("henn_RRR",b["henn"]["RRR"]),
            ("eeres_RRR",b["eeres"]["RRR"]),
            ("two_programme_pooled_RRR",b["exploratory_two_programme_mKH"]["pooled_RRR"])):
            if not math.isclose(float(a[key]),float(val),rel_tol=1e-11,abs_tol=1e-12):
                raise ValueError(f"STOP: source-derived output changed for {a['deleted_order']}: {key}")
        expected_ci=a["two_programme_mKH_t1_95"]
        observed_ci=b["exploratory_two_programme_mKH"]["mKH_t1_95"]
        if len(expected_ci)!=2 or len(observed_ci)!=2 or not all(
            math.isclose(float(x),float(y),rel_tol=1e-11,abs_tol=1e-12)
            for x,y in zip(expected_ci,observed_ci)):
            raise ValueError("STOP: meta uncertainty drift")
    return True


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("results/alan_drift_four_order_leaveoneout_20261124.json"))
    p.add_argument("--frozen",type=Path,default=None,help="Require all source-derived numerical outputs match this pinned snapshot")
    a=p.parse_args()
    output=run(base.source_original_xlsx())
    if a.frozen is not None:
        verify_frozen(output,json.loads(a.frozen.read_text(encoding="utf-8")))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("ORDER_CATCH_RATIOS",json.dumps({k:{o:round(v["lit_over_dark_catch_RR"],6) for o,v in t.items()}
        for k,t in output["per_order_catch_ratios_not_drift_density"].items()},ensure_ascii=False))
    print("FOUR_ORDER_SENSITIVITY",json.dumps([
       {"drop":d["deleted_order"],"Henn":d["henn"]["RRR"],"EERES":d["eeres"]["RRR"],
        "k2_RRR":d["exploratory_two_programme_mKH"]["pooled_RRR"],
        "k2_interval":d["exploratory_two_programme_mKH"]["mKH_t1_95"]}
       for d in output["comparisons"]],ensure_ascii=False))

if __name__=="__main__":
    main()
