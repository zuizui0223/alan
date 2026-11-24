#!/usr/bin/env python3
"""ALAN x crayfish x post-treatment week, original 16-flume 2022 experiment.

Ogbeide et al. 2026, Functional Ecology, DOI 10.1111/1365-2435.70335
Author-source DOI 10.5281/zenodo.14516568.

Input is the ORIGINAL-SOURCE-DERIVED flume-level aggregation of 192 author
emergence observations (two sampled weeks, six repeated trap-day records per
flume per week). The source RDS sha256 and R-to-CSV sha256 are pinned in
source_manifest.csv. This analysis replays the 32 derived flume-week rows,
NOT a fresh access to the 322MB author ZIP.

True treatment assignment unit: 16 physically independent stream-bank
mesocosms, 4 per 2x2 cell. Same 16 flumes revisited in W1 and W4.
Neither week is a true pre-intervention baseline, since light treatments were
already in place in both weeks. The six trap records are SUBSAMPLES.

Target post-hoc factorial contrast:
  [log(Emerge_L+cray/Emerge_cray) - log(Emerge_L/Emerge_control)]_W4
 -[log(Emerge_L+cray/Emerge_cray) - log(Emerge_L/Emerge_control)]_W1.

Paired flume within-cell (W4,W1) covariance is retained in the
multivariate delta SE. Only approximate Welch-Satterthwaite t df at very
small treatment n. Randomization of 16 full mesocosms is the assumption
behind a descriptive balanced-label permutation sharp-null sensitivity.
Do not claim a predation rate, predator-driven mechanism, plant fitness,
permanent abundance effect, causal mediation, or third pooled meta k.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import t

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"database/alan_experiments/ogbeide_sources/flume_emergence_rates.csv"
MANIFEST=ROOT/"database/alan_experiments/ogbeide_sources/source_manifest.csv"
WEEK4=ROOT/"database/alan_experiments/ogbeide_sources/week4_light_effects.csv"

EXPECTED_BLOB="a774508cdbe30641e97de61a7f55d437e15b9aec"
SOURCE_ORIGINAL_HASH={
    "emergence.rds":"289e9dc874724ca765aa2e8314f975ef7cf7b2cb688bdb62e8f1be9bc50d1776",
    "emergence_export.csv":"520c0ff129d68992b7e18f8c51ba54bec63bb13ed4c438f569c903c9e4c46d46",
}
CONDITIONS=("alan","control","cray","alan+cray")
WEEKS=("Week 1","Week 4")
N_PER_TREATMENT=4

def source_verify_and_read():
    original=SOURCE.read_bytes()
    blob=hashlib.sha1(f"blob {len(original)}\0".encode()+original).hexdigest()
    if blob!=EXPECTED_BLOB:
        raise ValueError("STOP: source-derived flume/week archive Git blob changed")
    with MANIFEST.open(encoding="utf-8",newline="") as fh:
        rows=list(csv.DictReader(fh))
    manifest={x["file"]:x["sha256"] for x in rows}
    if manifest!=SOURCE_ORIGINAL_HASH or len(rows)!=2:
        raise ValueError("STOP: original source RDS and export hash manifest changed")
    with SOURCE.open(encoding="utf-8",newline="") as fh:
        src=list(csv.DictReader(fh))
    if len(src)!=32:
        raise ValueError("STOP: source derived flume-week rows not 32")
    d={}
    for x in src:
        week=x["week"]; treatment=x["treat"]; flume=x["flume"]
        if week not in WEEKS or treatment not in CONDITIONS:
            raise ValueError("STOP: unexpected source treatment or week")
        rate=float(x["mean_rate"])
        capture=float(x["capture_sum"])
        if (x["trap_records"]!="6" or not flume.startswith("f")
                or not math.isfinite(rate) or rate<=0
                or not math.isclose(rate,capture/24,rel_tol=1e-10,abs_tol=1e-9)):
            raise ValueError("STOP: source trap exposure, positive rate or sample key mismatch")
        k=(flume,week)
        if k in d:raise ValueError("STOP: duplicated flume/week")
        d[k]=(treatment,rate)
    ids=sorted({f for f,w in d})
    if len(ids)!=16 or len(d)!=32:
        raise ValueError("STOP: not sixteen independent original mesocosms")
    arms={c:[] for c in CONDITIONS}
    for flume in ids:
        tr1,v1=d[flume,"Week 1"]
        tr4,v4=d[flume,"Week 4"]
        if tr1!=tr4:raise ValueError("STOP: original plot treatment changed between sampled weeks")
        arms[tr1].append((flume,v1,v4))
    if any(len(x)!=N_PER_TREATMENT for x in arms.values()):
        raise ValueError("STOP: four randomized mesocosms per 2x2 treatment required")
    effects={}
    with WEEK4.open(encoding="utf-8",newline="") as fh:
        checks={x["effect_id"]:x for x in csv.DictReader(fh)}
    if set(checks)!={"Ogbeide2026_week4_no_crayfish","Ogbeide2026_week4_with_crayfish"}:
        raise ValueError("STOP: original source Week4 paired effect ledger changed")
    for suffix,num,den in (("no_crayfish","alan","control"),
                           ("with_crayfish","alan+cray","cray")):
        observed=np.mean([r[2] for r in arms[num]])/np.mean([r[2] for r in arms[den]])
        reference=float(checks[f"Ogbeide2026_week4_{suffix}"]["ratio"])
        if not math.isclose(observed,reference,rel_tol=1e-11,abs_tol=1e-11):
            raise ValueError("STOP: true 4-flume Week4 comparison does not match original prior audit")
        effects[suffix]=reference
    return arms,blob,effects

def interaction_week(means,week):
    i=week
    return math.log(means["alan+cray"][i]/means["cray"][i])-math.log(means["alan"][i]/means["control"][i])

def data_summary(arms,permutation_count=60000,seed=20261124):
    means={c:np.mean(np.array([(row[1],row[2]) for row in arms[c]],float),axis=0)
           for c in CONDITIONS}
    week1=interaction_week(means,0);week4=interaction_week(means,1)
    change=week4-week1

    group_variance={}
    for cond in CONDITIONS:
        u=np.array([(x[1],x[2]) for x in arms[cond]],float)
        # Influence for the (log mean Week4 - log mean Week1) within
        # this physical treatment arm. Carries within-flume covariance.
        scores=(u[:,1]-means[cond][1])/means[cond][1] - (u[:,0]-means[cond][0])/means[cond][0]
        group_variance[cond]=float(scores.var(ddof=1)/4)
    v=sum(group_variance.values())
    se=math.sqrt(v)
    df=v*v/sum(g*g/3 for g in group_variance.values())
    q=float(t.ppf(.975,df))
    ci=[math.exp(change-q*se),math.exp(change+q*se)]
    p=float(2*t.sf(abs(change/se),df))

    # The source authors report random assignments to four balanced cells.
    # This is an approximate, conditional sharp-null permutation sensitivity,
    # not a preregistered design-based test; no blocking inferred.
    if permutation_count<1000 or permutation_count>300000:
        raise ValueError("STOP: bounded permutation count required")
    arrays=np.concatenate([
        np.array([(row[1],row[2]) for row in arms[cond]],float)
        for cond in CONDITIONS
    ],axis=0)
    rng=np.random.default_rng(seed)
    permutations=np.argsort(rng.random((permutation_count,16)),axis=1)
    permmeans=arrays[permutations].reshape(permutation_count,4,4,2).mean(axis=2)
    # axis 1 = (alan,control,cray,alan+cray)
    cross=np.log(permmeans[:,3,:]/permmeans[:,2,:])-np.log(permmeans[:,0,:]/permmeans[:,1,:])
    null_stat=cross[:,1]-cross[:,0]
    p_random=(1+int(np.count_nonzero(np.abs(null_stat)>=abs(change))))/(permutation_count+1)
    random_mc_se=math.sqrt(p_random*(1-p_random)/(permutation_count+1))

    return {
        "doi":"10.1111/1365-2435.70335",
        "data_doi":"10.5281/zenodo.14516568",
        "source_derived_file_git_blob":EXPECTED_BLOB,
        "n_true_independent_light_randomization_flumes":16,
        "n_flumes_per_2_by_2_cell":4,
        "measured_timepoints":list(WEEKS),
        "both_timepoints_after_light_on":True,
        "measured_endpoint":"Chironomidae emerged per m2 per recorded trap day, same 6 trap records per flume-week",
        "arm_mean_emergence_week1_week4":{k:list(map(float,v)) for k,v in means.items()},
        "no_crayfish_ALAN_vs_dark_RR_by_week":{
            "Week 1":float(means["alan"][0]/means["control"][0]),
            "Week 4":float(means["alan"][1]/means["control"][1])},
        "crayfish_present_ALAN_vs_dark_RR_by_week":{
            "Week 1":float(means["alan+cray"][0]/means["cray"][0]),
            "Week 4":float(means["alan+cray"][1]/means["cray"][1])},
        "light_by_crayfish_interaction_RRR":{
            "Week 1":math.exp(week1),
            "Week 4":math.exp(week4),
        },
        "interaction_on_original_log_scale":{
            "Week 1":week1,"Week 4":week4},
        "week4_minus_week1_change_in_ALAN_x_crayfish_log_interaction":change,
        "triple_interaction_ratio_of_RRRs":math.exp(change),
        "three_way_within_flume_delta_sampling_se":se,
        "four_cell_sampling_variance_contributions":group_variance,
        "welch_satterthwaite_df":df,
        "exploratory_paired_flume_t95_triple_ratio":ci,
        "exploratory_two_sided_t_p":p,
        "conditional_2_by_2_random_label_permutation_sharp_null":{
            "permutations_with_equal_cell_sizes":permutation_count,
            "seed":seed,
            "approx_two_sided_p":p_random,
            "Monte_Carlo_se":random_mc_se,
            "not_an_exact_pairwise_stratified_assignment_test":True,
            "potential_untested_spatial_blocking":True,
        },
        "status":"SOURCE_DERIVED_POSTHOC_2x2xTIME_DESCRIPTIVE_NOT_META",
        "independent_meta_k_increment":0,
        "limitations":[
            "Only chironomid emergence, not total aquatic adult assemblage or flying trap attraction",
            "The source had 16 true mesocosms assigned once; three emergence traps per flume and weeks were repeated subsamples",
            "Week1 and Week4 are both exposed and cannot estimate before-versus-after light introduction",
            "Week4 is not Week6 final experiment, source has no Week6 emergence",
            "No direct predator stomach contents/attack frequency tested here; crayfish presence is a manipulated treatment not a measured causal mediating predation path",
            "The permutation is a conditional sharp-null sensitivity under an assumed fully exchangeable 16-unit assignment; not replacement for original allocation protocol",
            "Posthoc outcome and week selection, weak statistical precision with four mesocosms per cell",
            "Study alone is one physical programme and cannot supply an additional independent meta effect from two times or 4 treatment groups",
            "Cannot be combined with EERES in-stream drift, flying trap catch or wet-plant reproductive fitness endpoints",
        ],
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("results/alan_ogbeide2026_2x2xweek_emergence_20261124.json"))
    parser.add_argument("--permutations",type=int,default=60000)
    args=parser.parse_args()
    arms,sha,checks=source_verify_and_read()
    report=data_summary(arms,permutation_count=args.permutations)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("OGBEIDE_SOURCE_REPLAY",json.dumps({
        "author_flume_rates_blob_verified":sha,
        "source_week4_rr_check":checks,
        "true_n":report["n_true_independent_light_randomization_flumes"],
        "interaction_week1":report["light_by_crayfish_interaction_RRR"]["Week 1"],
        "interaction_week4":report["light_by_crayfish_interaction_RRR"]["Week 4"],
        "triple_RRR":report["triple_interaction_ratio_of_RRRs"],
        "t95":report["exploratory_paired_flume_t95_triple_ratio"],
        "t_p":report["exploratory_two_sided_t_p"],
        "permutation":report["conditional_2_by_2_random_label_permutation_sharp_null"],
        "meta_k_added":0
    },ensure_ascii=False))

if __name__=="__main__":
    main()
