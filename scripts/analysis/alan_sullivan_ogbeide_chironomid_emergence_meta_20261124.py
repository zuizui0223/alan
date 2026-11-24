#!/usr/bin/env python3
"""Source-verified, prefit-defined exploratory 2-programme Chironomidae emergence.

Programme S: Sullivan2013/2019, Ohio, three habitat blocks, four artificial
light configurations each contrasted to one shared dark within each block.
Programme O: Ogbeide2022/2026, Germany, four independent ALAN and four
independent dark stream-bank flumes with no crayfish, two after-light sampling
weeks. Source and location are distinct from Sullivan.

The estimand is within-programme log-lit/dark of author-labelled Chironomidae
caught in aquatic emergence traps. Similar stage/taxon; light colour,
intensity, trapping method, dates and habitat are NOT transportable.
Not a formally registered systematic/meta or causal mediator.
"""
import argparse
import importlib.util
import json
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SOURCE_SCRIPT=ROOT/"scripts/analysis/alan_sullivan2019_chironomid_source_gate_20261124.py"
OGBEIDE_SCRIPT=ROOT/"scripts/analysis/alan_ogbeide2026_emergence_light_crayfish_time_20261124.py"
T2=4.302652729911275
T1=12.706204736432095

def import_named(name,filename):
    spec=importlib.util.spec_from_file_location(name,filename)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

sullivan=import_named("original_sullivan2019_family_emergence",SOURCE_SCRIPT)
ogbeide=import_named("original_ogbeide2026_family_emergence",OGBEIDE_SCRIPT)

def independent_habitat_block_estimates(source):
    if source.get("status")!="SOURCE_LABEL_RETRIEVED_NO_META_EFFECT":
        raise ValueError("STOP: no author-original Sullivan source label gate")
    if source["original_2013_author_TSV_sha256"]!=sullivan.DIGEST:
        raise ValueError("STOP: original Sullivan source hash drift")
    if source["chironomid_source_labels"]!=["Chironomidae"]:
        raise ValueError("STOP: Chironomidae not a single explicitly tagged taxon label")
    if source["original_stations"]!=15 or source["independent_habitat_blocks"]!=3:
        raise ValueError("STOP: original source physical block count drift")
    expected={"DRK":356,"LEDI":1162,"LEDH":1295,"HPSI":1457,"HPSH":690}
    if source["chironomid_count_by_treatment"]!=expected:
        raise ValueError("STOP: original author Chironomidae total or arm changed")
    group={}
    for r in source["source_block_treatment_counts"]:
        g=str(r["group"]);arm=r["light"]
        if g not in ("1","2","3") or arm not in sullivan.ARMS or (g,arm) in group:
            raise ValueError("STOP: source block/light duplicated or missing")
        if int(r["four_original_dates"])!=4:
            raise ValueError("STOP: missing matched source four sample dates")
        total=int(r["all_aquatic_capture"])
        chiro=int(r["explicit_chironomid_recorded_capture"])
        if not 0<=chiro<=total:
            raise ValueError("STOP: invalid family share in original emergent source counts")
        group[g,arm]=chiro
    if len(group)!=15:
        raise ValueError("STOP: original wetland block x light identities missing")
    # Per-block geometric average of four dependent, same-dark station,
    # original light menus. Shared dark is not counted as four experiments.
    def selection(lights):
        blocklog=[]
        for g in ("1","2","3"):
            d=group[g,"DRK"]
            if d<=0:raise ValueError("STOP: actual dark Chiro group zero; no pseudocount")
            vals=[]
            for tr in lights:
                n=group[g,tr]
                if n<=0:raise ValueError("STOP: actual lit Chiro group zero; no pseudocount")
                vals.append(math.log(n/d))
            blocklog.append(statistics.mean(vals))
        mu=statistics.mean(blocklog)
        se=statistics.stdev(blocklog)/math.sqrt(3)
        return {"light_menu":list(lights),"three_block_lnRR":blocklog,
                "true_block_n":3,"logRR":mu,"SE":se,"RR":math.exp(mu),
                "descriptive_t2_95":[math.exp(mu-T2*se),math.exp(mu+T2*se)]}
    all_lights=("LEDI","LEDH","HPSI","HPSH")
    original=selection(all_lights)
    variants={
        "all_four_equal_colour":original,
        "LED_only_equal_colour":selection(("LEDI","LEDH")),
        "HPS_only_equal_colour":selection(("HPSI","HPSH")),
    }
    for lamp in all_lights:
        variants["only_"+lamp]=selection((lamp,))
    return original,variants,{f"{g}/{t}":group[g,t] for g in ("1","2","3") for t in sullivan.ARMS}

def independent_flume_estimates(arms):
    if set(arms)!={"alan","control","cray","alan+cray"} or any(len(x)!=4 for x in arms.values()):
        raise ValueError("STOP: original 16 flumes 2x2 four per cell changed")
    if len({r[0] for arm in arms.values() for r in arm})!=16:
        raise ValueError("STOP: original flume identity collision")
    # Original source no-crayfish eight true flumes; same flume Week1/Week4.
    # Mean within-flume two weeks, then mean across independent flumes.
    def one_week_choice(kind):
        if kind not in ("both_equal_weeks","Week 1","Week 4"):
            raise ValueError("STOP: unknown prespecified time contrast")
        treatment={}
        for c in ("alan","control"):
            u=[]
            for _,w1,w4 in arms[c]:
                value=statistics.mean((w1,w4)) if kind=="both_equal_weeks" else (w1 if kind=="Week 1" else w4)
                if not math.isfinite(value) or value<=0:
                    raise ValueError("STOP: selected true-flume positive rate missing")
                u.append(float(value))
            treatment[c]=u
        lit,unlit=treatment["alan"],treatment["control"]
        ml,md=statistics.mean(lit),statistics.mean(unlit)
        logrr=math.log(ml/md)
        v_lit=statistics.variance(lit)/(4*ml*ml)
        v_dark=statistics.variance(unlit)/(4*md*md)
        variance=v_lit+v_dark
        return {
            "original_flume_n_lit":4,"original_flume_n_dark":4,
            "sampled_week_choice":kind,
            "source_lit_mean":ml,"source_dark_mean":md,
            "flume_unit_sampling_variance_delta":variance,
            "logRR":logrr,"SE":math.sqrt(variance),"RR":math.exp(logrr),
            "flume_mean_values_by_treatment":treatment,
            "not_independent_week_replicates":True,
        }
    variants={week:one_week_choice(week) for week in ("both_equal_weeks","Week 1","Week 4")}
    primary=variants["both_equal_weeks"]
    prior_week4=ogbeide.source_verify_and_read()[2]["no_crayfish"]
    if not math.isclose(variants["Week 4"]["RR"],prior_week4,rel_tol=1e-12,abs_tol=1e-12):
        raise ValueError("STOP: prior true-flume Week4 light-control effect not reproduced")
    return primary,variants

def k2_meta(s,o):
    y=[float(s["logRR"]),float(o["logRR"])]
    v=[float(s["SE"])**2,float(o["SE"])**2]
    if min(v)<=0 or not all(math.isfinite(z) for z in y+v):
        raise ValueError("STOP: missing programme-level sampling variance")
    tau2=max(0,((y[0]-y[1])**2-v[0]-v[1])/2)
    w=[1/(z+tau2) for z in v]
    mu=sum(z*q for z,q in zip(y,w))/sum(w)
    q=sum(ww*(yy-mu)**2 for ww,yy in zip(w,y))
    se=math.sqrt(max(1.0,q)/sum(w))
    return {
        "k_independent_programmes":2,
        "REML_tau2_closed_form_k2":tau2,
        "mKH_residual_scale_floor":max(1,q),
        "pooled_logRR":mu,"pooled_RR":math.exp(mu),
        "pooled_logRR_SE_mKH":se,
        "modified_HK_t1_95":[math.exp(mu-T1*se),math.exp(mu+T1*se)],
        "df_for_meta_interval":1,
        "not_confirmatory":True,
    }


def independent_programme_difference_descriptive(s,o):
    """Contrast separate studies, not a causal wetland × ALAN moderator.

    The wetland estimate is the geometric mean across n=3 paired source
    blocks, and the no-crayfish flume estimate is log of arithmetic mean
    ratio across n=4+4 physical assignments. These are not identical
    statistical targets, so a t2 CI is a descriptive uncertainty alarm,
    not a confirmatory test of effect modification.
    """
    y=float(s["logRR"])-float(o["logRR"])
    se=math.hypot(float(s["SE"]),float(o["SE"]))
    if not (se>0 and math.isfinite(y)):
        raise ValueError("STOP: invalid independent-programme comparison")
    t2=4.302652729911275
    tval=abs(y/se)
    p_t2=1-tval/math.sqrt(tval*tval+2)
    return {
        "source_two_programme_RR_of_RRs":math.exp(y),
        "source_lnRR_difference":y,
        "uncertainty_SE_from_two_independent_source_programmes":se,
        "descriptive_conservative_t2_95_RR_of_RRs":[math.exp(y-t2*se),math.exp(y+t2*se)],
        "descriptive_t2_two_sided_p":p_t2,
        "fraction_naive_sampling_variance_from_Sullivan_n3":s["SE"]**2/se**2,
        "not_causal_wetland_vs_flume_effect_modifier":True,
        "not_confirmatory_cross_study_difference":True,
        "method_warning":"Different source estimands, independent sources but n3 wetland blocks and n8 flume plots; t2 interval is an intentionally conservative descriptive reference, not a general meta-regression effect."
    }


def compute(source_original,original_flume_arms):
    wet,wetvariants,units=independent_habitat_block_estimates(source_original)
    flume,flumevariants=independent_flume_estimates(original_flume_arms)
    primary=k2_meta(wet,flume)
    # Predeclared leave-one-physical-unit-out diagnostics are point-only.
    # With 2 wetland blocks remaining, do NOT attach population CIs.
    block_logs=wet["three_block_lnRR"]
    wet_loo={
        f"group_{i+1}":math.exp(statistics.mean(
            y for j,y in enumerate(block_logs) if j!=i))
        for i in range(3)
    }
    # Remove one real independently assigned Ogbeide flume, retaining
    # within-flume Week1/Week4 and unchanged no-crayfish design.
    original_lit=flume["flume_mean_values_by_treatment"]["alan"]
    original_dark=flume["flume_mean_values_by_treatment"]["control"]
    flume_loo={}
    for treatment,vals in (("alan",original_lit),("control",original_dark)):
        for i in range(4):
            if treatment=="alan":
                r=statistics.mean(v for j,v in enumerate(vals) if j!=i)/statistics.mean(original_dark)
            else:
                r=statistics.mean(original_lit)/statistics.mean(v for j,v in enumerate(vals) if j!=i)
            flume_loo[f"{treatment}_index_{i+1}"]=r
    diagnostics={
        "Sullivan_n3_group_lnRR":block_logs,
        "Sullivan_n3_group_ratios":[math.exp(y) for y in block_logs],
        "Sullivan_groups_positive_count":sum(y>0 for y in block_logs),
        "Sullivan_leave_one_group_out_RR_point_only":wet_loo,
        "Ogbeide_eight_flume_leave_one_out_RR_point_only":flume_loo,
        "Ogbeide_Week1_RR":flumevariants["Week 1"]["RR"],
        "Ogbeide_Week4_RR":flumevariants["Week 4"]["RR"],
        "between_programme_ratio_of_LitDark_RR_point_only":wet["RR"]/flume["RR"],
        "note":"Post hoc diagnostics share original physical source programmes, no additional k or causal moderator"
    }
    sensitivity={
        f"wetland_{sk}_Ogbeide_all_weeks":k2_meta(sv,flume)
        for sk,sv in wetvariants.items()

    }
    sensitivity.update({
        f"wetland_all_four_Ogbeide_{k}":k2_meta(wet,val)
        for k,val in flumevariants.items()
    })
    return {
        "status":"SOURCE_AUTHENTICATED_EXPLORATORY_CHIRONOMIDAE_EMERGENCE_TWO_PROGRAMME_K2",
        "design_protocol":"docs/ALAN_CHIRONOMID_EMERGENCE_TWO_PROGRAMME_PROTOCOL_20261124.md",
        "source_study_dois":["10.1002/eap.1821","10.1111/1365-2435.70335"],
        "source_original_sha256":[sullivan.DIGEST,ogbeide.SOURCE_ORIGINAL_HASH["emergence.rds"]],
        "shared_endpoint_description":"Source-labelled Chironomidae emergent adult catches/rates; wetland natural field versus artificial freshwater flume",
        "Sullivan_source_units_three_group_by_treatment_chiro":units,
        "Sullivan_color_menu_frozen_primary":wet,
        "Sullivan_source_colour_sensitivities":wetvariants,
        "Ogbeide_no_crayfish_eight_flumes_frozen_primary":flume,
        "Ogbeide_week_sensitivities":flumevariants,
        "exploratory_k2_random_effects_meta":primary,
        "two_programme_differential_effect_size_descriptive":independent_programme_difference_descriptive(wet,flume),
        "leave_one_physical_unit_out_point_sensitivities":diagnostics,
        "contrasts_if_changing_source_colour_or_sampling_week":sensitivity,
        "new_study_programmes_discovered_today":0,
        "new_exploratory_endpoint_meta_lane":1,
        "added_to_existing_108_effect_catalogue":False,
        "changed_existing_four_meta_effects":False,
        "is_final_PRISMA":False,
        "caveats":[
            "Sullivan samples 4 colour/dose arms that share each block's single dark station, Ogbeide neutral white LED only",
            "Source trap gear/area and effort differ; within-study log ratio cancels common effort only under equality of original source sampling dates/units",
            "Original Sullivan n=3 paired wetland blocks is too small for reliable population inference",
            "Ogbeide n=4 lit + 4 unlit only in no-crayfish cells, original wetland and flume environments differ",
            "Both studies have original source record and taxon labels verified, but human independent raw rereading not completed",
            "The two-study random effects variance and t1 interval are highly unstable; no typical universal ALAN emergence effect established",
            "No model of wet versus dry, plant reproduction, predation pathway or causal mediator; week sensitivity is post-treatment",
            "Taxonomic identity same family, not necessarily species/functional guild; site/gear and ALAN dose modifiers not estimable k2",
        ],
    }


def verify_frozen(actual, frozen):
    """Fail-closed source-based replay of the literal k2 paper-facing values."""
    if frozen.get("status")!="FROZEN_AFTER_ORIGINAL_SOURCE_REPLAY_EXPLORATORY_NOT_FINAL_SYSTEMATIC":
        raise ValueError("STOP: no prior source-verified frozen snapshot")
    if actual.get("status")!="SOURCE_AUTHENTICATED_EXPLORATORY_CHIRONOMIDAE_EMERGENCE_TWO_PROGRAMME_K2":
        raise ValueError("STOP: actual result not source authenticated")
    pairs=[
        (actual["Sullivan_color_menu_frozen_primary"]["RR"], frozen["programme_estimates"][0]["RR"]),
        (actual["Sullivan_color_menu_frozen_primary"]["SE"], frozen["programme_estimates"][0]["SE_logRR"]),
        (actual["Ogbeide_no_crayfish_eight_flumes_frozen_primary"]["RR"], frozen["programme_estimates"][1]["RR"]),
        (actual["Ogbeide_no_crayfish_eight_flumes_frozen_primary"]["SE"], frozen["programme_estimates"][1]["SE_logRR"]),
        (actual["exploratory_k2_random_effects_meta"]["REML_tau2_closed_form_k2"], frozen["exploratory_meta"]["REML_tau2"]),
        (actual["exploratory_k2_random_effects_meta"]["pooled_RR"], frozen["exploratory_meta"]["pooled_RR"]),
        (actual["exploratory_k2_random_effects_meta"]["pooled_logRR_SE_mKH"], frozen["exploratory_meta"]["SE_logRR_mKH"]),
    ]
    pairs.extend(zip(actual["exploratory_k2_random_effects_meta"]["modified_HK_t1_95"],frozen["exploratory_meta"]["t1_95"]))
    contrast=actual.get("two_programme_differential_effect_size_descriptive")
    saved=frozen.get("two_programme_difference_descriptive")
    if not contrast or not saved:
        raise ValueError("STOP: missing cross-programme effect-size difference sensitivity")
    for key in (
        "source_two_programme_RR_of_RRs",
        "source_lnRR_difference",
        "uncertainty_SE_from_two_independent_source_programmes",
        "descriptive_t2_two_sided_p",
    ):
        pairs.append((contrast[key],saved[key]))
    pairs.extend(zip(contrast["descriptive_conservative_t2_95_RR_of_RRs"],
                     saved["descriptive_conservative_t2_95_RR_of_RRs"]))
    if (not contrast["not_causal_wetland_vs_flume_effect_modifier"] or
        not contrast["not_confirmatory_cross_study_difference"]):
        raise ValueError("STOP: two studies cannot identify causal wetland moderator")

    if any(not math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-12) for a,b in pairs):
        raise ValueError("STOP: published source-based family-level emergence estimate or interval drift")
    if frozen.get("Sullivan2013_original_Dryad_TSV_SHA256")!=sullivan.DIGEST:
        raise ValueError("STOP: Sullivan raw bytes identity drift")
    if frozen.get("Ogbeide2022_original_Zenodo_RDS_SHA256")!=ogbeide.SOURCE_ORIGINAL_HASH["emergence.rds"]:
        raise ValueError("STOP: Ogbeide original RDS identity drift")
    if frozen.get("Ogbeide2022_derived_source_flume_rates_git_blob")!=ogbeide.EXPECTED_BLOB:
        raise ValueError("STOP: Ogbeide flume aggregate Git blob identity drift")
    if actual["exploratory_k2_random_effects_meta"]["k_independent_programmes"]!=2:
        raise ValueError("STOP: pseudo programme replication attempted")
    prior=frozen.get("within_programme_physical_unit_diagnostics")
    diag=actual.get("leave_one_physical_unit_out_point_sensitivities")
    if not prior or not diag:
        raise ValueError("STOP: missing physical-unit leave-one-out source audit")
    indexed=[
        (diag["Sullivan_n3_group_ratios"],prior["Sullivan_3_group_RRs"]),
        ([diag["Sullivan_leave_one_group_out_RR_point_only"]["group_"+str(i)] for i in (1,2,3)],
         prior["Sullivan_leave_one_group_out_RRs"]),
        ([diag["Ogbeide_eight_flume_leave_one_out_RR_point_only"][f"{arm}_index_{i}"]
         for arm in ("alan","control") for i in (1,2,3,4)],
         prior["Ogbeide_8_flume_leave_one_out_RRs"]),
    ]
    for observed,pinned in indexed:
        if len(observed)!=len(pinned) or not all(math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12)
                                                   for a,b in zip(observed,pinned)):
            raise ValueError("STOP: physical-unit robustness/source effect drift")
    for key,reference in (
        ("Ogbeide_Week1_RR",prior["Ogbeide_Week1_RR"]),
        ("Ogbeide_Week4_RR",prior["Ogbeide_Week4_RR"]),
        ("between_programme_ratio_of_LitDark_RR_point_only",prior["between_programme_ratio_of_RRs_point"]),
    ):
        if not math.isclose(float(diag[key]),float(reference),rel_tol=1e-10,abs_tol=1e-12):
            raise ValueError("STOP: source conditional-week or cross-programme point drift")
    return True


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("results/alan_sullivan_ogbeide_chironomid_emergence_k2_20261124.json"))
    p.add_argument("--frozen",type=Path,default=None,help="Require byte-authenticated source results to reproduce pinned paper-facing estimates")
    args=p.parse_args()
    source=sullivan.analyze(sullivan.get_original())
    flumes,flumeblob,_=ogbeide.source_verify_and_read()
    report=compute(source,flumes)
    if args.frozen is not None:
        verify_frozen(report,json.loads(args.frozen.read_text(encoding="utf-8")))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("TWO_PROGRAMME_CHIRONOMID_EMERGENCE",json.dumps({
        "sourceS_original_sha_verified":True,
        "sourceO_flume_blob_verified":flumeblob,
        "Sullivan_all_four":{k:report["Sullivan_color_menu_frozen_primary"][k] for k in ("RR","SE","descriptive_t2_95")},
        "Sullivan_by_lamp":{name:{"RR":val["RR"],"SE":val["SE"]}
             for name,val in report["Sullivan_source_colour_sensitivities"].items()},
        "Ogbeide_no_cray_all_weeks":{k:report["Ogbeide_no_crayfish_eight_flumes_frozen_primary"][k] for k in ("RR","SE")},
        "pooled_k2":report["exploratory_k2_random_effects_meta"],
        "opposing_programme_ratio_comparison":report["two_programme_differential_effect_size_descriptive"],
        "physical_unit_loo":report["leave_one_physical_unit_out_point_sensitivities"],
        "sensitivity":{k:v["pooled_RR"] for k,v in report["contrasts_if_changing_source_colour_or_sampling_week"].items()},
        "existing_four_meta_effects_unchanged":True
    },ensure_ascii=False))

if __name__=="__main__":
    main()
