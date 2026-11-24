#!/usr/bin/env python3
"""Exploratory k=2 meta of source-matched four insect orders in aquatic drift.

Canonical 108-effect database and existing stage meta results remain untouched.
Per-source effects are ln ratio-of-ratios using the ratio of study arm means
across original independent spatial units; *not* mean within-stream log odds.
Sampling variance is a paired, four-group multivariate delta approximation.

For two studies, closed-form intercept-only REML tau² follows from the
restricted likelihood of y1-y2 ~ N(0, s1²+s2²+2*tau²). Modified Knapp-Hartung
uses max(1, Q/(k-1)) to prevent a variance interval narrower than
the ordinary random-effects interval. df=k-1=1; do not use normal 1.96.
"""
import argparse
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
INPUT=HERE/"two_programme_four_order_source_estimates_20261124.json"
T_CRIT_DF1=12.706204736432095
T_CRIT_DF4=2.7764451051977987
T_CRIT_DF5=2.570581835636305

def review_source_lock(meta,replay):
    if replay.get("status")!="SOURCE_AUTHENTICATED_EXPLORATORY_TWO_PROGRAMME_COMPARISON_NOT_ACCEPTED_META":
        raise ValueError("STOP: no original-source replay")
    if replay.get("meta_k_added")!=0 or not replay.get("mkh_or_reml_pooled_not_fitted"):
        raise ValueError("STOP: unexpected source replay state")
    values={ "Henn2012Texas_five_stream_pairs":replay["henn"],
             "EERES2021Germany_six_header_tanks":replay["eeres"] }
    for r in meta["studies"]:
        a=values.get(r["programme_id"])
        if a is None:
            raise ValueError("STOP: programme lost in replay")
        for frozen,fresh in (("ln_rrr","four_order_log_composition_RRR"),
                             ("se_paired_physical_unit_delta","delta_SE_shared_physical_units"),
                             ("rrr","four_order_composition_RRR")):
            if not math.isclose(float(r[frozen]),float(a[fresh]),rel_tol=1e-11,abs_tol=1e-12):
                raise ValueError(f"STOP: author source changed {r['programme_id']} {frozen}")
        if r["n_physical_units"]!=a["n_independent_physical_units"]:
            raise ValueError("STOP: physical independence count changed")
    if replay["eeres"].get("source_sha256")!=meta["eeres_dryad_workbook_sha256"]:
        raise ValueError("STOP: source hash mismatch")
    return True


def synthesize(meta):
    studies=meta.get("studies")
    if len(studies)!=2:
        raise ValueError("STOP: exactly 2 independent source programmes are available")
    ids=[r["programme_id"] for r in studies]
    if set(ids)!={"Henn2012Texas_five_stream_pairs","EERES2021Germany_six_header_tanks"}:
        raise ValueError("STOP: changed or double-counted programme identity")
    if meta.get("numerator_orders")!=["Ephemeroptera","Trichoptera"] or meta.get("denominator_orders")!=["Diptera","Coleoptera"]:
        raise ValueError("STOP: estimand taxa not identical")
    if meta.get("canonical_108_effect_catalogue_changed") or meta.get("existing_three_stage_meta_k_changed"):
        raise ValueError("STOP: original results must be preserved")
    if meta.get("status")!="TWO_SOURCE_AUTHENTICATED_PROGRAMME_ESTIMATES_FOR_EXPLORATORY_SYNTHESIS_NOT_FINAL_PRISMA":
        raise ValueError("STOP: no trusted source-extraction status")
    y=[float(r["ln_rrr"]) for r in studies]
    s=[float(r["se_paired_physical_unit_delta"]) for r in studies]
    n=[int(r["n_physical_units"]) for r in studies]
    if n != [5,6] or any(not math.isfinite(x) or x<=0 for x in s):
        raise ValueError("STOP: source granularity or SE invalid")
    if any(not math.isclose(math.exp(e),float(r["rrr"]),rel_tol=1e-11) for e,r in zip(y,studies)):
        raise ValueError("STOP: study-scale RR and lnRR disagree")
    v=[z*z for z in s]

    # With k=2, the difference is the only residual contrast, and the exact
    # restricted log likelihood is equivalent to a single normal variance.
    diff=y[0]-y[1]
    vdiff=v[0]+v[1]
    tau2=max(0.0,(diff**2-vdiff)/2.0)
    w=[1/(z+tau2) for z in v]
    mu=sum(y_i*w_i for y_i,w_i in zip(y,w))/sum(w)
    q=sum(w_i*(y_i-mu)**2 for y_i,w_i in zip(y,w))  # df1
    q_mkh=max(1.0,q)
    se_mu=math.sqrt(q_mkh/sum(w))
    lo,hi=math.exp(mu-T_CRIT_DF1*se_mu),math.exp(mu+T_CRIT_DF1*se_mu)

    # Henn author's family sum is 2 more than the printed Total for the
    # San Marcos light arm. Bound point estimate if all 2 were either
    # E+T or D+C. Leave source SE unchanged in this *point-only* sensitivity.
    henn=studies[0]
    if not math.isclose(y[0],math.log((737*688)/(402*1113)),rel_tol=1e-12):
        raise ValueError("STOP: unexpected Henn source groups")
    bounds=[]
    for A_lit,B_lit in ((735,402),(737,402),(737,400)):
        altered=math.log((A_lit*688)/(B_lit*1113))
        y_alt=[altered,y[1]]
        tau_alt=max(0,((y_alt[0]-y_alt[1])**2-vdiff)/2)
        w_alt=[1/(z+tau_alt) for z in v]
        mu_alt=sum(z*q for z,q in zip(y_alt,w_alt))/sum(w_alt)
        bounds.append(math.exp(mu_alt))

    return {
        "status":"EXPLORATORY_TWO_PROGRAMME_META_NOT_SYSTEMATIC_CONFIRMATORY",
        "effect_definition":meta["estimand"],
        "studies":[{"programme":ids[i],"physical_n":n[i],"lnRRR":y[i],"SE":s[i],"RRR":math.exp(y[i])} for i in range(2)],
        "k_independent_programmes":2,
        "method":"intercept-only normal random effects, REML tau2 closed for k=2; modified Knapp-Hartung floor=1, t(df=1)",
        "pooled_lnRRR":mu,"pooled_RRR":math.exp(mu),
        "pooled_mkh_95_interval":[lo,hi],
        "pooled_lnRRR_se_mkh":se_mu,
        "t_reference_df":1,
        "reml_tau2":tau2,
        "reml_tau":math.sqrt(tau2),
        "Q_unadjusted_df1":(diff**2)/vdiff,
        "Q_reml_df1":q,
        "mkh_variance_factor":q_mkh,
        "two_count_Henn_reallocation_pooled_point_range":[min(bounds),max(bounds)],
        "input_but_not_added_to_main_effect_catalogue":True,
        "new_stage_meta_distinct_from_existing_three":True,
        "limitations":[
            "k=2 gives very unstable between-study variance and very broad t1 interval; cannot infer a universal positive ecological effect",
            "Henn authors' Table 1 family row totals disagree with printed Total by 2 in San Marcos light; prose totals also differ",
            "Henn source taxonomy manually transcribed into six orders, not a byte-exact archived data file; EERES genuine original XLSX SHA checked",
            "Henn uses one lit and one dark night per natural stream with high-intensity incandescent light; EERES experimental flumes use repeated weeks and different LED exposure",
            "Study-level paired physical-block delta sampling variances are approximations and have tiny n=5/6",
            "Taxa restricted to E+T vs Diptera+Coleoptera; not all drifting arthropods, population abundances, predator diet or wetland plants",
            "Observational effects on wet-preferring plant populations and pollen limitation are not causally explained",
            "Author reports and partial systematic screen not yet fully independently adjudicated; this is not a finalized PRISMA systematic meta-analysis",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--frozen",type=Path,default=INPUT)
    ap.add_argument("--original-replay",type=Path,default=None,
                    help="Optional newly downloaded EERES source replay file for strict source-lock check")
    ap.add_argument("--output",type=Path,default=Path("results/alan_four_order_two_programme_meta_20261124.json"))
    args=ap.parse_args()
    frozen=json.loads(args.frozen.read_text(encoding="utf-8"))
    if args.original_replay is not None:
        original=json.loads(args.original_replay.read_text(encoding="utf-8"))
        review_source_lock(frozen,original)
    outcome=synthesize(frozen)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(outcome,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("TWO_PROGRAMME_META",json.dumps({
        "status":outcome["status"],
        "k":outcome["k_independent_programmes"],
        "RRR":outcome["pooled_RRR"],
        "95_mKH":outcome["pooled_mkh_95_interval"],
        "tau2":outcome["reml_tau2"],
        "source_replay_checked":args.original_replay is not None,
    },ensure_ascii=False))


if __name__=="__main__":
    main()
