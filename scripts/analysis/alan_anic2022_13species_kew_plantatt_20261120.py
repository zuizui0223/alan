#!/usr/bin/env python3
"""Original Anic 2022 13 meadow plants x independent UK PLANTATT Ellenberg F.

This is a NEW trait-source sensitivity on existing original 665 plant-trough
biomass measurements already authenticated and summarized on this research
branch (original 54 independent physical troughs, source 2022 Figshare PDF).
The new independently maintained RBG Kew 2022 BIFloraExplorer raw CSV
provides PLANTATT (Hill et al. 2004, corrected 2008) indicator F for 12
exact binomials + one Agrostis tenuis -> A. capillaris validated by Kew
POWO. NEVER use fuzzy name matching or replace unknown traits with dry.

Source BIFloraExplorer Git blob hash (SHA-1 Git object) is pinned.
Only plant-level log biomass ratio across L1/L2/L3 original *unmapped*
treatment codes is analyzed. Need source author code key to identify dark,
amber, and white LED. Absent biomass source cells are not filled with zeros.

Null: species labels for F shuffled across all 12 complete-case plant taxa,
calculate maximum absolute Spearman rho across three planned light-code
comparisons. 4,999 permutations and leave-physical-plant-one-out diagnostics.
Plant taxa within one physical 54-trough experiment are not independently
randomized ALAN treatments; no climate/hydrology treatment is present.
The one high-wet F=8 taxon means MOIST>7 contrast remains unsupported.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,random
from pathlib import Path
from scipy.stats import spearmanr

SCHEMA="alan_anic2022_13species_PLANTATT_F_original_lightcode_source_v1"
RAW_URL="https://raw.githubusercontent.com/RBGKew/BIFloraExplorer/master/data-raw/BI_MAIN.csv"
RAW_GIT_BLOB="e6a1da929af81f370d2a1b79177b489e3e573a06"
SOURCE_JSON="results/all_species_forest/biotic_interactions/anic2022_original_13species_eiv_moisture_lightcode_20261120.json"
SOURCE_PDF_SHA="baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0"
K_GIT_SOURCE="https://github.com/RBGKew/BIFloraExplorer/blob/master/data-raw/BI_MAIN.csv"
K_POWO_SYN="https://powo.science.kew.org/taxon/386415-1"
SYN={"Agrostis tenuis":"Agrostis capillaris"}
KEYS=("L2/L1","L3/L1","L3/L2")
SEED=20261120
BOOT=4999
def require(x,msg):
    if not x:raise RuntimeError("ORIGINAL_Kew_PLANTATT_SOURCE_STOP: "+msg)

def git_blob_hash(raw):
    return hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()

def pull_kew_binary():
    import urllib.request
    req=urllib.request.Request(RAW_URL,headers={
      "User-Agent":"ALAN original Kew source ecological audit"})
    with urllib.request.urlopen(req,timeout=80) as f:
        raw=f.read(3<<20)
    require(git_blob_hash(raw)==RAW_GIT_BLOB,
            "Kew original BIFloraExplorer raw CSV git object sha mismatch")
    return raw

def trait_lookup(raw,min_rows=3000):
    reader=csv.DictReader(io.StringIO(raw.decode("latin1")))
    require({"taxon_name_binom","F_PLANTATT","taxonomic_status",
             "aggregate","F_Doring"}.issubset(reader.fieldnames or []),
            "Kew raw plant naming/F columns changed")
    m={}
    count=0
    for r in reader:
        count+=1
        name=r["taxon_name_binom"].strip()
        f=r["F_PLANTATT"].strip()
        if f in ("","NA","NaN","nan","x","X"):continue
        try:val=float(f)
        except ValueError:continue
        require(math.isfinite(val) and 1<=val<=12,
                "Kew original F out of accepted Ellenberg range")
        if name in m:
            require(m[name]["F_PLANTATT"]==val,
                    "Kew same binomial conflicting F values")
        m[name]={"F_PLANTATT":val,
                  "Kew_taxonomic_status":r["taxonomic_status"].strip(),
                  "Kew_aggregate_label":r["aggregate"].strip(),
                  "F_Doring":r["F_Doring"].strip()}
    require(count>=min_rows,"Kew source expected full flora rows")
    return m,count

def original_source(path):
    z=json.loads(Path(path).read_text())
    require(z["schema"]==
            "alan_Anic2022_13plant_original_EIV_exact_source_verified_compact_v1",
            "original author 13 source result schema drift")
    require(z["original_source_pdf_sha256"]==SOURCE_PDF_SHA,
            "Anic 2022 author-original PDF SHA changed")
    require(z["original_plant_species"]==13 and
            z["source_treatment_physical_mesocosms"]==54,
            "original 13 plants 54 random lighting troughs changed")
    obs=z["observed_13_species_source_light_code_logratios"]
    require(len(obs)==13 and len({x["species"] for x in obs})==13,
            "original Anic plant source observations differ")
    return obs

def rho(x,y):
    if len(x)<4 or len(set(x))<2 or len(set(y))<2:return None
    out=float(spearmanr(x,y).statistic)
    return out if math.isfinite(out) else None

def describe(obs,traits):
    out=[]
    for r in obs:
        orig=r["species"]
        verified_name=SYN.get(orig,orig)
        info=traits.get(verified_name)
        out.append({
           "source_name":orig,
           "original_biomass_source_pdf_SHA256":SOURCE_PDF_SHA,
           "independent_Kew_name":verified_name if info else None,
           "name_match_method":"authoritative_Kew_POWO_synonym" if orig in SYN
                else ("original_exact_binomial" if info else "no_match"),
           "F_PLANTATT":info["F_PLANTATT"] if info else None,
           "F_Doring":info["F_Doring"] if info else None,
           "Kew_taxon_aggregate_flag":info["Kew_aggregate_label"] if info else None,
           "source_lightcode_log_biomass_ratio":r["codes"],
           "lightcode_spectral_key_not_identified":True,
        })
    return out

def infer(rows,reps=BOOT):
    allF=[x for x in rows if x["F_PLANTATT"] is not None]
    eligible=[x for x in allF if all(x["source_lightcode_log_biomass_ratio"].get(c) is not None
                                    for c in KEYS)]
    wet=sum(x["F_PLANTATT"]>7 for x in allF)
    low=sum(x["F_PLANTATT"]<=7 for x in allF)
    results={
      "original_plant_names":len(rows),"source_F_nonmissing":len(allF),
      "source_all_three_code_complete_plant_names":len(eligible),
      "moisture_above_7_original_plant_names":wet,
      "moisture_at_or_below_7_plant_names":low,
      "source_F_value_counts":{str(f):sum(x["F_PLANTATT"]==f for x in allF)
                                for f in sorted({x["F_PLANTATT"] for x in allF})},
      "source_F_PLANTATT_range":[min(x["F_PLANTATT"] for x in allF),
                                 max(x["F_PLANTATT"] for x in allF)],
      "binary_wet_moisture_test_supported":wet>=4 and low>=4,
      "original_source_figshare_lightcode_to_color_mapping_verified":False,
      "source_phylogeny_or_plot_covariance_modeled":False,
      "same_source_54_light_trophic_troughs":True,
      "status":"EXPLORATORY_ORDINAL_MOISTURE_GRADIENT_ONLY"
    }
    require(len(eligible)>=10,"source comparable all-three light contrasts too sparse")
    f=[x["F_PLANTATT"] for x in eligible]
    effect={k:[x["source_lightcode_log_biomass_ratio"][k] for x in eligible] for k in KEYS}
    point={k:rho(f,y) for k,y in effect.items()}
    require(all(x is not None for x in point.values()),"source Spearman undefined")
    # Pre-specified entire-species F permutation across three contrasts.
    rng=random.Random(SEED)
    maxnull=[]
    singles={k:[] for k in KEYS}
    for _ in range(reps):
        pf=f.copy();rng.shuffle(pf)
        values={k:rho(pf,y) for k,y in effect.items()}
        require(all(x is not None for x in values.values()),
                "source null Spearman undefined")
        maxnull.append(max(abs(x) for x in values.values()))
        for k,val in values.items():singles[k].append(abs(val))
    perm_adj={k:(1+sum(x>=abs(r)-1e-12 for x in maxnull))/(reps+1)
              for k,r in point.items()}
    perm_nom={k:(1+sum(x>=abs(r)-1e-12 for x in singles[k]))/(reps+1)
              for k,r in point.items()}
    # Species identity and light-code sensitivity; high-wet one-species
    # and POWO synonym/aggregate exclusions must be transparent.
    loo={}
    for k,y in effect.items():
        loo[k]={
          x["source_name"]:rho([v for i,v in enumerate(f) if i!=j],
                              [v for i,v in enumerate(y) if i!=j])
          for j,x in enumerate(eligible)
        }
    without_high={}
    nohigh=[x for x in eligible if x["F_PLANTATT"]<=7]
    for k in KEYS:
        without_high[k]=rho([x["F_PLANTATT"] for x in nohigh],
                              [x["source_lightcode_log_biomass_ratio"][k] for x in nohigh])
    no_syn=[x for x in eligible if x["name_match_method"]=="original_exact_binomial"]
    no_agg=[x for x in eligible if x["Kew_taxon_aggregate_flag"] not in ("agg","in agg")]
    def getrho(sub):
        return {k:rho([x["F_PLANTATT"] for x in sub],
                      [x["source_lightcode_log_biomass_ratio"][k] for x in sub])
                 for k in KEYS}
    results.update({
      "source_spearman_rho":point,
      "source_species_label_permutation_nominal_p":perm_nom,
      "max_three_preplanned_comparisons_adjusted_p":perm_adj,
      "permutations":reps,
      "source_plant_leave_one_out_rho":loo,
      "without_single_high_wet_Lotus_pedunculatus_source_rho":without_high,
      "without_single_high_wet_source_species_n":len(nohigh),
      "exclude_authoritative_POWO_synonym_source_rho":getrho(no_syn),
      "exclude_agg_source_rho":getrho(no_agg),
      "no_realized_habitat_water_or_reproductive_endpoint":True,
      "no_ALAN_causal_moisture_moderator_identified":True,
    })
    return results

def run(kew_raw,source_path,reps=BOOT):
    traits,count=trait_lookup(kew_raw)
    original=original_source(source_path)
    rows=describe(original,traits)
    results=infer(rows,reps)
    return {
      "schema":SCHEMA,
      "source_Kew_GitBlob_SHA1":RAW_GIT_BLOB,
      "source_Kew_nomenclature_url":K_GIT_SOURCE,
      "source_synonym_validation_url":K_POWO_SYN,
      "source_original_Kew_BIFloraExplorer_rows":count,
      "source_original_data_Anic_PDF_SHA256":SOURCE_PDF_SHA,
      "original_plant_species_F_and_biomass_ratios":rows,
      "inference_and_support":results,
      "not_source_isolated_independently_illuminated_plants":True,
      "not_a_meta_analysis":True,
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source",default=SOURCE_JSON)
    a.add_argument("--out",required=True)
    args=a.parse_args()
    r=run(pull_kew_binary(),args.source)
    p=Path(args.out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(r,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("ANIC13_KEW_PLANTATT_SOURCE_TRAIT_COVERAGE_AND_MAX3_NULL "+json.dumps(
      {"source_rows":r["source_original_Kew_BIFloraExplorer_rows"],
       "traits":r["original_plant_species_F_and_biomass_ratios"],
       "result":r["inference_and_support"]},sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
