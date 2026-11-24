#!/usr/bin/env python3
"""Audit covariance-unidentified Lu JPE2026 night-vs-open reproductive contrast.

Use primary publisher-S4b source-derived estimates already extracted in ALAN:
DOI 10.1093/jpe/rtag160. Original S4b contrast direction Unlit minus Lit.

B = [Lit - Unlit]_(Open) - [Lit - Unlit]_(Nocturnal)
  = [Unlit - Lit]_(Nocturnal) - [Unlit - Lit]_(Open).

This B is on the UNKNOWN REPORTED MODEL SCALE. It is never a seed/ovule
percentage-point change, lnRR, or causal compensation fraction. Unknown
shared-block covariance rho in [-1, 1] induces a sharp marginal SE envelope.
A 1.96*z display is a sensitivity heuristic, NOT a valid 95% confidence
interval in the small-number-of-light-installations field design.
"""
import argparse
import csv
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parents[2]
SOURCE=HERE/"database/alan_experiments/lu_jpe_sources/table_S4b.csv"
EXPECTED_DOI="10.1093/jpe/rtag160"
SOURCE_ORIGINAL_CONTRAST="Unlit - Lit"

# Source-value guards for intended two seed comparisons.
ORIGINAL_NUMERIC={
    ("P. kansuensis","Seed set","Open"):(0.0445,0.0319,92.8),
    ("P. kansuensis","Seed set","Diurnal"):(-0.0686,0.0401,96.8),
    ("P. kansuensis","Seed set","Nocturnal"):(0.1332,0.0429,95.5),
    ("E. regelii","Seed set","Open"):(-0.0094,0.0308,80.9),
    ("E. regelii","Seed set","Diurnal"):(0.0062,0.0333,81.5),
    ("E. regelii","Seed set","Nocturnal"):(-0.0040,0.0340,80.9),
}
RESPONSE_TYPES=("Seed set","Pollen number")
ALLOWED_ACCESS=("Open","Diurnal","Nocturnal")

def read_source():
    with SOURCE.open(newline="",encoding="utf-8") as f:
        raw=list(csv.DictReader(f))
    if len(raw)!=12:raise ValueError("STOP: primary original S4b is not twelve contrasts")
    species=response=None
    table={}
    for r in raw:
        species=r["Species"].strip() or species
        response=r["Response"].strip() or response
        access=r["Pollination treatment"].strip()
        if species not in {"P. kansuensis","E. regelii"} or response not in RESPONSE_TYPES or access not in ALLOWED_ACCESS:
            raise ValueError("STOP: unexpected original S4b experimental cell")
        if r["Contrast"].strip()!=SOURCE_ORIGINAL_CONTRAST or r["doi"]!=EXPECTED_DOI:
            raise ValueError("STOP: reported effect orientation or DOI changed")
        key=(species,response,access)
        if key in table:raise ValueError("STOP: duplicate reported contrast")
        effect=float(r["Estimate"])
        se=float(r["SE"])
        df=float(r["df"])
        if not all(map(math.isfinite,(effect,se,df))) or se<=0 or df<=0:
            raise ValueError("STOP: missing original model uncertainty")
        table[key]={"effect_unlit_minus_lit":effect,"model_SE":se,"model_df":df,"reported_p":float(r["p value"])}
    if len(table)!=12:raise ValueError("STOP: incomplete factorial design")
    for k,(v,se,df) in ORIGINAL_NUMERIC.items():
        r=table[k]
        if not all(math.isclose(float(r[n]),y,rel_tol=0,abs_tol=1e-9)
                   for n,y in zip(("effect_unlit_minus_lit","model_SE","model_df"),(v,se,df))):
            raise ValueError(f"STOP: original seed treatment contrast changed {k}")
    return table

def contrast_envelope(table,species,response):
    o=table[species,response,"Open"]
    n=table[species,response,"Nocturnal"]
    # Source Unlit-Lit orientation: night minus open makes the intended B.
    B=n["effect_unlit_minus_lit"]-o["effect_unlit_minus_lit"]
    so,sn=o["model_SE"],n["model_SE"]
    min_se=abs(sn-so)
    max_se=sn+so
    z_critical=1.96
    se_rho_zero=math.hypot(sn,so)
    threshold_rho=None
    if so*sn>0:
        threshold_rho=(so*so+sn*sn-(abs(B)/z_critical)**2)/(2*so*sn)
    display={
        "species":species,
        "outcome":response,
        "B_definition":"(Lit-Unlit)_Open - (Lit-Unlit)_Nocturnal, source model contrast units",
        "B_point_on_original_reported_model_scale":B,
        "night_unlit_minus_lit":n["effect_unlit_minus_lit"],
        "open_unlit_minus_lit":o["effect_unlit_minus_lit"],
        "SE_night_from_report":sn,
        "SE_open_from_report":so,
        "unknown_rho_interval":[-1,1],
        "SE_bound_cauchy_schwarz":[min_se,max_se],
        "hypothetical_rho_zero_SE":se_rho_zero,
        "worst_case_normal_z_1p96_envelope":[B-z_critical*max_se,B+z_critical*max_se],
        "rho_needed_for_naive_nominal_z1p96_nonzero":threshold_rho,
        "source_model_scale_unverified":True,
        "physical_light_plot_covariance_unknown":True,
        "proper_design_based_95CI_estimable":False,
        "absolute_seed_fraction_or_pp_change":None,
        "adopt_into_meta":False
    }
    return display

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=Path("results/alan_lu2026_access_covariance_envelope_20261124.json"))
    args=ap.parse_args()
    table=read_source()
    records=[contrast_envelope(table,spp,response)
             for spp in ("P. kansuensis","E. regelii")
             for response in ("Seed set","Pollen number")]
    data={
        "study":"Lu et al. 2026 JPE",
        "doi":EXPECTED_DOI,
        "input_original_report_table":"database/alan_experiments/lu_jpe_sources/table_S4b.csv",
        "source_recovered_supplement":"publisher DOCX SHA256 c9cc7ea119c6a99bcfee2b6956c76988c73cf776f20944f5570d00101c3e3abe",
        "status":"SOURCE_REPORTED_COEFFICIENTS_ONLY_UNKNOWN_JOINT_COVARIANCE",
        "study_experimental_programmes_added":0,
        "already_known_in_ALAN":True,
        "source_12_rows_not_12_independent_effects":True,
        "estimates":records,
        "key_inference":"For P kansuensis, comparing open versus nocturnal ALAN effects is algebraically possible but SE of difference cannot be identified without the 5 paired field-block joint covariance and source scale. Do not claim daytime compensation from point estimate or separate p values.",
        "independent_light_installations_verified":False,
        "light_by_access_omnibus_in_source":"P kansuensis seed S4a reported F=6.11, P=.0031: does not test specific open vs nocturnal B directly.",
        "method_warning":"Cauchy-Schwarz bounds on missing covariance, 1.96 envelope only normal-approximation sensitivity; not small-sample confirmatory CI.",
        "main_108_catalogue_unchanged":True,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":data["status"],
                      "P_kansuensis_seed_set":records[0],
                      "E_regelii_seed_set":records[2]},ensure_ascii=False))

if __name__=="__main__":
    main()
