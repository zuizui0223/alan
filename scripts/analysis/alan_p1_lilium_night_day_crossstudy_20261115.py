#!/usr/bin/env python3
"""P1 H4: genuine Kendall/Nicholson source rows for Izu vs Honshu Lilium.

Compare *night-only/day-only* seed:ovule ratios, NOT the primary day/open
marginal-night-loss contrast: Morinaga 2009's third arm is archived as
hand cross-pollination with an ambiguous original figure caption.
Only 2 papers, 4 study-year conditions. No ALAN exposure measured here.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics
from pathlib import Path

SHA="1d9c4c778206584130308b99ba21d902b0ccf1047118dfcf52abaec9d652e2aa"
DOIS={"Nakajima_2018":"10.1111/1442-1984.12222",
      "Morinaga_2009":"10.1007/s10144-008-0097-1"}
ARMS={"day pollination","night pollination"}
EXPECTED={
 ("Morinaga_2009",2003):(0.418103,0.318966,13,19),
 ("Morinaga_2009",2004):(0.573276,0.508621,3,8),
 ("Morinaga_2009",2005):(0.517241,0.633621,14,13),
 ("Nakajima_2018",2016):(0.03,0.49,4,9)
}
def must(condition,reason):
    if not condition:raise ValueError("P1_LILIUM_SOURCE_GATE: "+reason)
def valid_float(v):
    f=float(v)
    must(math.isfinite(f) and f>=0,"nonnegative finite value")
    return f
def extract(raw):
    must(hashlib.sha256(raw).hexdigest()==SHA,"original extraction-table SHA drift")
    rows=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    must(len(rows)==1652,"source table row count drift")
    pool=[x for x in rows if x["study_ID"] in DOIS]
    must(len(pool)==15,"source Nakajima+Morinaga paper row coverage drift")
    grouped={}
    for r in pool:
        st=r["study_ID"];y=int(r["year_start"]);k=(st,y)
        must(r["DOI"]==DOIS[st],"original article DOI changed")
        must(r["treatment_effectiveness_metric"] in ("seed set","fruit set"),
             "other outcome added")
        if r["treatment_effectiveness_metric"]!="seed set":continue
        grouped.setdefault(k,{})[r["treatment_condition"]]=r
    must(set(grouped)==set(EXPECTED),"source study/year rows changed")
    out=[]
    for (st,y),a in sorted(grouped.items()):
        must(ARMS.issubset(a),"original day+night access missing")
        # Third arm is NOT safely interpretable as open for Morinaga2009.
        if st=="Morinaga_2009":
            must(set(a)==ARMS|{"hand pollination"},"mainland third arm no longer hand-pollination")
            must("hand cross" in a["hand pollination"]["treatment_comments"],
                 "original hand/open conflict lost; revisit Fig3 and methods")
        else:must(set(a)==ARMS|{"open pollination"},"island open arm changed")
        d=a["day pollination"];n=a["night pollination"]
        must(d["effectiveness_units"] in ("seed/ovule","seed:ovule")
             and n["effectiveness_units"] in ("seed/ovule","seed:ovule"),
             "shared seed/ovule outcome invalid")
        D=valid_float(d["effectiveness_value"]);N=valid_float(n["effectiveness_value"])
        sD=valid_float(d["effect_error"]);sN=valid_float(n["effect_error"])
        nd=int(float(d["sample_size"]));nn=int(float(n["sample_size"]))
        must(d["effect_error_units"]=="SD" and n["effect_error_units"]=="SD",
             "source error was not sample SD")
        must(D>0 and N>0 and nd>0 and nn>0,"log contrast requires positive means/n")
        exp=EXPECTED[(st,y)]
        must(abs(D-exp[0])<1e-9 and abs(N-exp[1])<1e-9 and (nd,nn)==exp[2:],
             "original year-specific figure-extracted mean or n changed")
        lr=math.log(N/D)
        # Exploratory delta-method SE; ignores figure extraction error,
        # within-paper year correlations and original randomization method.
        variance=sN*sN/(nn*N*N)+sD*sD/(nd*D*D)
        se=math.sqrt(variance)
        out.append(dict(study_ID=st,publication_DOI=DOIS[st],year=y,
                        island_status="Izu-Oshima" if st=="Nakajima_2018" else "Honshu",
                        daytime_seed_ovule=D,nighttime_seed_ovule=N,
                        day_sample_n=nd,night_sample_n=nn,day_sd=sD,night_sd=sN,
                        log_night_day=lr,night_day_ratio=N/D,
                        approximate_independent_arm_sampling_se=se,
                        original_source_extraction=d["extraction_type"],
                        note="means and SDs extracted from published figure, not raw plant observations"))
    return out

def summarize(items):
    island=[r for r in items if r["study_ID"]=="Nakajima_2018"]
    mainland=[r for r in items if r["study_ID"]=="Morinaga_2009"]
    must(len(island)==1 and len(mainland)==3,"paper/year independence gate")
    isl=island[0];i_log=isl["log_night_day"]
    yearly=[r["log_night_day"] for r in mainland]
    m_mean=statistics.mean(yearly)
    delta=i_log-m_mean
    m_var=sum(r["approximate_independent_arm_sampling_se"]**2 for r in mainland)/9
    approx_se=math.sqrt(isl["approximate_independent_arm_sampling_se"]**2+m_var)
    iv_w=[1/r["approximate_independent_arm_sampling_se"]**2 for r in mainland]
    iv_log=sum(w*r["log_night_day"] for w,r in zip(iv_w,mainland))/sum(iv_w)
    deltas={str(r["year"]):math.exp(i_log-statistics.mean(
            j["log_night_day"] for j in mainland if j["year"]!=r["year"]))
            for r in mainland}
    # Explicit *rounding-only* sensitivity. Digital figure values are not known
    # accurate to ±0.005: this is a minimal numerical perturbation check.
    low=(isl["nighttime_seed_ovule"]-.005)/(isl["daytime_seed_ovule"]+.005)/math.exp(m_mean)
    high=(isl["nighttime_seed_ovule"]+.005)/(isl["daytime_seed_ovule"]-.005)/math.exp(m_mean)
    return {
      "endpoint":"night-only / day-only proportion of seeds per ovule",
      "design":"two distinct papers, island n=1 study-year, mainland n=3 study-years within ONE paper",
      "independent_publications":2,"independent_author_publication_clusters":2,
      "years_are_not_three_independent_mainland_studies":True,
      "island_night_day_ratio":isl["night_day_ratio"],
      "mainland_year_log_effects":{str(r["year"]):r["log_night_day"] for r in mainland},
      "mainland_equal_year_mean_log_night_day":m_mean,
      "mainland_geometric_average_night_day_ratio":math.exp(m_mean),
      "island_minus_mainland_log_night_day_equal_year":delta,
      "night_day_ratio_of_ratios_island_over_mainland":math.exp(delta),
      "sampling_only_approximate_log_difference_se":approx_se,
      "sampling_only_approximate_ratio95_WALD_not_full_CI":[
        math.exp(delta-1.96*approx_se),math.exp(delta+1.96*approx_se)],
      "iv_weighted_mainland_year_sensitivity_ratio_of_ratios":math.exp(i_log-iv_log),
      "leave_one_mainland_year_out_ratio_of_ratios":deltas,
      "figure_rounding_plus_minus_0p005_island_only_ratio_of_ratios":[low,high],
      "island_day_vs_night_access_hours":9,
      "mainland_day_vs_night_access_hours":12,
      "island_mesh_exclusion_study_durations_days":2,
      "mainland_mesh_exclusion_full_blooming_period":True,
      "same_study_island_mainland_paired":False,
      "same_year_island_mainland_paired":False,
      "source_mainland_control_really_open":"UNVERIFIED_CONFLICTING_FIGURE_AND_METHODS",
      "night_dependence_from_day_over_open_estimable_in_mainland":False,
      "ALAN_exposure_measured":False,
      "causal_island_adaptation_established":False,
      "source_figure_extraction_error_quantified":False,
      "source_arm_independence_assumed_for_wald":True,
      "inferential_warning":"Delta-method CI is a fragile within-figure sampling-only illustration; not a cross-study causal island-effect CI or a valid two-publication meta-analysis."
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    source=args.source.read_bytes()
    rows=extract(source)
    result={"schema":"ALAN_P1_Lilium_Izu_Honshu_2003to2016_cross_study_figures_v1",
            "source_sha256":SHA,"source_extraction_table_rows":1652,
            "year_arms":rows,"contrast":summarize(rows)}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    must(not args.out.exists(),"refuse overwrite")
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_P1_LILIUM_REAL_SOURCE",json.dumps({"rows":len(rows),
        "geographic_log_effect":result["contrast"]["island_minus_mainland_log_night_day_equal_year"],
        "ratio_of_ratios":result["contrast"]["night_day_ratio_of_ratios_island_over_mainland"],
        "samples_only_interval":result["contrast"]["sampling_only_approximate_ratio95_WALD_not_full_CI"],
        "mainland_control":result["contrast"]["source_mainland_control_really_open"],
        "ALAN_measured":False},ensure_ascii=False))
if __name__=="__main__":main()
