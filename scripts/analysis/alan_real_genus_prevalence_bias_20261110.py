#!/usr/bin/env python3
"""Source-locked genus and documentation-prevalence audit of ALAN × plant records.

Read only original 93,209,351 photo-ID source-derived accepted plant×GSHHG
island×0.05° site×year rows, and MD5-audited 1992–2024 prior-year PCNL.
The GloBI moth label means *published visit/pollinates record*, NOT effective
nocturnal moth dependence. Missing records are unknown, not biological absence.
Never infer cultivation status or regional native status from photo metadata.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr
from scipy.stats import t

SITE_SHA = "440a8313d71168e0a4dfd2f41a7bef45623abb068b07e3282026f1c35dccbd2c"
GLOBI_SHA = "a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994"
LIGHT = "pcnl_1992_to_previous_year_mean_log1p"
KEY = ["island", "latbin", "lonbin", "year"]
SCHEMA = "ALAN_real93m_genus_and_source_ubiquity_bias_v1"


def require(ok, why):
    if not ok:
        raise ValueError("ALAN_GENUS_SOURCE_BIAS: " + str(why))


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_sites(path):
    require(sha256(path) == SITE_SHA, "Original signed 33-year PCNL input changed")
    d = pd.read_csv(path)
    require(len(d) == 85318 and not d.duplicated(KEY).any(),
            "Original 85,318 photo-site/year denominators changed")
    d["island"] = d.island.astype(str)
    d = d.loc[d.historic_ntl_complete.astype(str).str.lower().eq("true") &
              d.pcnl_years_observed.eq(d.pcnl_years_expected) &
              d.pcnl_years_expected.eq(d.year-1992) &
              d.all_accepted_recorded_species.ge(3) &
              d.all_accepted_observation_uuids.ge(5)].copy()
    q = d.groupby(["island", "year"], observed=True)[LIGHT].agg(["size", "min", "max"])
    keep = q.loc[q["size"].ge(2) & (q["max"]-q["min"]).ge(.02)].index.get_level_values(0).unique()
    d = d.loc[d.island.isin(keep) & d[LIGHT].notna()].copy()
    require(len(d) == 77647 and d.island.nunique() == 562,
            "Original within-island ALAN contrast support changed")
    d["log_obs"] = np.log1p(d.all_accepted_observation_uuids)
    return d


def load_labels(globi):
    require(sha256(globi) == GLOBI_SHA, "Frozen source GloBI study checksum changed")
    g = pd.read_csv(globi, compression="gzip", dtype=str).fillna("")
    require(len(g) == 61967 and set(g.guild) ==
            {"moth_family_visit", "butterfly_family_visit"},
            "61,967 moth/butterfly GloBI source claims missing")
    g["moth"] = g.guild.eq("moth_family_visit")
    g["butterfly"] = g.guild.eq("butterfly_family_visit")
    g["moth_external"] = g.moth & ~g.study_title.str.contains(
        "iNaturalist", case=False, regex=False)
    f = g.groupby("plant_name", observed=True)[
        ["moth", "butterfly", "moth_external"]].any()
    studies = g.loc[g.moth_external & g.study_title.ne("")].groupby(
        "plant_name", observed=True).study_title.nunique()
    f["n_external_study_titles"] = studies.reindex(f.index, fill_value=0).astype(int)
    return f


def original_accepted_species(groups, sites, labels):
    frames = []
    for i in range(16):
        p = Path(groups) / ("allflora-island-group-%02d.sqlite" % i)
        receipt = p.with_suffix(".json")
        require(p.is_file() and receipt.is_file(), "Missing original source group %d" % i)
        a = json.loads(receipt.read_text())
        require(a.get("group") == i and
                a.get("schema") ==
                "ALAN_full93m_genuine_WCVP_allflora_moth_evidence_island_source_group_v1",
                "Wrong source group provenance")
        with sqlite3.connect(p) as con:
            x = pd.read_sql_query(
                "SELECT island,latbin,lonbin,year,species,"
                "moth_positive_photo_ids,butterfly_positive_photo_ids FROM flora", con)
        x["island"] = x.island.astype(str)
        x = x.merge(sites[KEY], on=KEY, how="inner", validate="many_to_one")
        frames.append(x)
    d = pd.concat(frames, ignore_index=True)
    d["moth_original"] = d.moth_positive_photo_ids.gt(0)
    d["butterfly_original"] = d.butterfly_positive_photo_ids.gt(0)
    duplicate = d.duplicated(KEY + ["species"], keep=False)
    if duplicate.any():
        agreements = d.loc[duplicate].groupby(KEY + ["species"], observed=True)[
            ["moth_original", "butterfly_original"]].nunique()
        require(agreements.to_numpy().max() <= 1,
                "Original source groups conflict on species GloBI evidence")
    d = d.drop_duplicates(KEY + ["species"]).copy()
    for k in ("moth", "butterfly", "moth_external", "n_external_study_titles"):
        d[k] = d.species.map(labels[k])
    for k in ("moth", "butterfly", "moth_external"):
        d[k] = d[k].fillna(False).astype(bool)
    d["n_external_study_titles"] = d.n_external_study_titles.fillna(0).astype(int)
    require(d.moth.eq(d.moth_original).all() and
            d.butterfly.eq(d.butterfly_original).all(),
            "Independent GloBI name join disagrees with 93m accepted source flags")
    d["genus"] = d.species.str.extract(r"^([A-Z][a-zA-Z-]+)\b", expand=False).fillna(
        "UNKNOWN_GENUS")
    counts = d.groupby(KEY, observed=True).size().reset_index(name="species_count")
    verify = sites[KEY + ["all_accepted_recorded_species"]].merge(
        counts, on=KEY, how="left", validate="one_to_one")
    require(len(d) == 2012925 and verify.species_count.notna().all() and
            verify.species_count.eq(verify.all_accepted_recorded_species).all(),
            "Source taxon×site records do not reconstruct 77,647 denominators")
    return d


def fixed_effects_regression(y, light, logobs, weights, island, fsets):
    """Weighted sparse FWL with true island-cluster covariance, 95% t CI."""
    y = np.asarray(y, dtype=float)
    light = np.asarray(light, dtype=float)
    logobs = np.asarray(logobs, dtype=float)
    w = np.asarray(weights, dtype=float)
    require(len(w) == len(y) and np.isfinite(w).all() and (w > 0).all(),
            "Invalid positive FWL weights")
    sw = np.sqrt(w)
    group_codes = [pd.factorize(f, sort=False)[0] for f in fsets]
    nc = sum(int(k.max()+1) for k in group_codes)
    idx = np.tile(np.arange(len(y), dtype=np.int32), len(group_codes))
    cols = []
    vals = []
    shift = 0
    for k in group_codes:
        norm = np.bincount(k, weights=w)
        cols.append(k + shift)
        vals.append(sw / np.sqrt(norm[k]))
        shift += len(norm)
    z = coo_matrix((np.concatenate(vals), (idx, np.concatenate(cols))),
                   shape=(len(y), nc)).tocsr()
    residualized = []
    iterations = []
    for series in (y, light, logobs):
        target = sw * series
        solve = lsmr(z, target, atol=1.e-8, btol=1.e-8, maxiter=1200)
        require(solve[1] in (1, 2), "Sparse source FE solver not converged")
        residualized.append((target-z @ solve[0])/sw)
        iterations.append(int(solve[2]))
    yy, xx, xo = residualized
    X = np.column_stack((xx, xo))
    xtwx = X.T @ (w[:, None] * X)
    require(np.linalg.matrix_rank(xtwx) == 2, "Light effect rank deficient")
    bread = np.linalg.inv(xtwx)
    beta = bread @ (X.T @ (w * yy))
    resid = yy-X @ beta
    codes = pd.factorize(island, sort=False)[0]
    n_islands = int(codes.max()+1)
    score = np.zeros((n_islands, 2))
    np.add.at(score, codes, X * (w * resid)[:, None])
    dof = len(y)-2-z.shape[1]
    require(n_islands >= 30 and dof > 0, "Insufficient island clusters or degrees")
    vcov = (n_islands/(n_islands-1)) * ((len(y)-1)/dof) * bread @ (score.T @ score) @ bread
    se = float(np.sqrt(max(float(vcov[0, 0]), 0)))
    q = float(t.ppf(.975, n_islands-1))
    slope = float(beta[0])
    return {"beta_pp": 100*slope,
            "95CI_pp": [100*(slope-q*se), 100*(slope+q*se)],
            "nominal_p": float(2*t.sf(abs(slope/se), n_islands-1)) if se > 0 else None,
            "n_site_years": None, "n_islands": n_islands,
            "FE_solver_iterations": iterations}


def species_model(species, sites, *, response="moth", genus_fe=True,
                  weight="island_equal"):
    counts = species.groupby(KEY, observed=True).size().reset_index(name="site_n")
    selected = sites.merge(counts, on=KEY, how="inner", validate="one_to_one")
    require(len(selected) >= 100 and selected.island.nunique() >= 30,
            "Species subset lacks enough geographic source sites")
    sub = species.merge(selected[KEY + ["site_n", LIGHT, "log_obs"]],
                        on=KEY, how="inner", validate="many_to_one")
    n_island = selected.groupby("island", observed=True).size()
    w = 1.0/sub.site_n.to_numpy(float)
    if weight == "island_equal":
        w = w/sub.island.map(n_island).to_numpy(float)
    elif weight != "site_equal":
        raise ValueError("Unsupported weighting")
    iy = sub.island.astype(str).str.cat(sub.year.astype(str), sep="|").to_numpy()
    fes = [iy] + ([sub.genus.to_numpy()] if genus_fe else [])
    fit = fixed_effects_regression(sub[response].astype(float),
           sub[LIGHT].to_numpy(float), sub.log_obs.to_numpy(float),
           w, sub.island.to_numpy(), fes)
    fit.update({"n_species_site_years": len(sub), "n_site_years": len(selected),
                "n_observed_genera": int(sub.genus.nunique()),
                "genus_FE": genus_fe, "island_x_year_FE": True, "weight": weight})
    return fit


def site_outcome_model(sites, count):
    x = sites.copy()
    x["response"] = 100.0 * np.asarray(count, dtype=float) / x.all_accepted_recorded_species
    # Source response is in percentage points; convert below to 0–1 before fit.
    w = 1.0/x.island.map(x.groupby("island", observed=True).size()).to_numpy(float)
    iy = x.island.astype(str).str.cat(x.year.astype(str), sep="|").to_numpy()
    v = fixed_effects_regression(x.response/100.0, x[LIGHT], x.log_obs,
                                 w, x.island.to_numpy(), [iy])
    v["n_site_years"] = len(x)
    return v


def run(sites_path, globi, groups, out):
    sites = load_sites(sites_path)
    flags = load_labels(globi)
    data = original_accepted_species(groups, sites, flags)
    taxa = data.drop_duplicates("species")[
        ["species", "genus", "moth", "moth_external", "butterfly", "n_external_study_titles"]]
    prevalence = data.groupby("species", observed=True).size().sort_values(
        ascending=False).rename("site_years")
    species_table = prevalence.reset_index().merge(taxa, on="species", validate="one_to_one")
    top1_n = int(np.ceil(len(species_table)*.01))
    top5_n = int(np.ceil(len(species_table)*.05))
    top1 = set(species_table.head(top1_n).species)
    top5 = set(species_table.head(top5_n).species)
    genus = taxa.groupby("genus", observed=True).agg(
        n_species=("species", "size"), n_moth=("moth", "sum"),
        n_butterfly=("butterfly", "sum")).reset_index()
    mixed = set(genus.loc[genus.n_species.ge(5) &
                genus.n_moth.ge(2) & (genus.n_species-genus.n_moth).ge(2), "genus"])
    data["moth_minus_butterfly"] = data.moth.astype(int)-data.butterfly.astype(int)
    ms = data.loc[data.genus.isin(mixed)]
    models = {}
    specs = [
       ("moth_no_genus_FE", data, "moth", False, "island_equal"),
       ("moth_genus_FE", data, "moth", True, "island_equal"),
       ("moth_non_iNat_genus_FE", data, "moth_external", True, "island_equal"),
       ("butterfly_genus_FE", data, "butterfly", True, "island_equal"),
       ("moth_minus_butterfly_paired_genus_FE", data, "moth_minus_butterfly", True, "island_equal"),
       ("mixed_180_genera_moth_genus_FE", ms, "moth", True, "island_equal"),
       ("mixed_180_genera_moth_non_iNat_genus_FE", ms, "moth_external", True, "island_equal"),
       ("moth_genus_FE_drop_top1pct_species", data.loc[~data.species.isin(top1)],
        "moth", True, "island_equal"),
       ("moth_genus_FE_drop_top5pct_species", data.loc[~data.species.isin(top5)],
        "moth", True, "island_equal"),
       ("moth_genus_FE_site_equal", data, "moth", True, "site_equal"),
    ]
    for name, frame, col, genus_fe, weight in specs:
        print("ALAN_GENUS_MODEL_RUNNING", name, len(frame), flush=True)
        models[name] = species_model(frame, sites, response=col,
                                      genus_fe=genus_fe, weight=weight)
    require(abs(models["moth_minus_butterfly_paired_genus_FE"]["beta_pp"] -
                (models["moth_genus_FE"]["beta_pp"]-
                 models["butterfly_genus_FE"]["beta_pp"])) < 1.e-7,
            "Moth versus butterfly must use identical source species site rows")
    require(abs(models["moth_no_genus_FE"]["beta_pp"] - .45722827881259) < .00002,
            "Source y and island×year FE did not reproduce verified base model")
    # Exact coefficient decomposition: same 562 islands / 77,647 cells and
    # denominator, with source top-1% species selected only by prevalence.
    group = data.assign(top1_positive=data.moth & data.species.isin(top1),
                        rest_positive=data.moth & ~data.species.isin(top1))
    agg = group.groupby(KEY, observed=True).agg(
        original=("moth", "sum"), top1=("top1_positive", "sum"),
        rest=("rest_positive", "sum")).reset_index()
    site = sites.merge(agg, on=KEY, validate="one_to_one")
    require(site.original.eq(site.moth_visitor_positive_recorded_species).all()
            and site.top1.add(site.rest).eq(site.original).all(),
            "Common-plant numerator does not conserve original moth response")
    contribution = {k: site_outcome_model(site, site[k]) for k in
                    ("original", "top1", "rest")}
    require(abs(contribution["original"]["beta_pp"] -
                contribution["top1"]["beta_pp"] -
                contribution["rest"]["beta_pp"]) < 1.e-8,
            "The documented GloBI source decomposition must be exact")
    require(abs(contribution["top1"]["beta_pp"] - .6318209378927077) < .00002,
            "Top prevalent taxa contribution not identical to original source")
    b = [0, 1, 4, 19, 99, 499, 10**9]
    names = ["1", "2-4", "5-19", "20-99", "100-499", "500+"]
    species_table["source_prevalence_bin"] = pd.cut(
        species_table.site_years, bins=b, labels=names)
    bins = species_table.groupby("source_prevalence_bin", observed=True).agg(
        n_species=("species", "size"), n_moth_positive=("moth", "sum"),
        n_external_moth=("moth_external", "sum")).reset_index()
    bins["moth_record_positive_pct"] = 100*bins.n_moth_positive/bins.n_species
    out = Path(out)
    require(not out.exists(), "Refuse overwrite of original empirical audit")
    out.mkdir(parents=True)
    species_table.to_csv(out/"accepted_species_source_ubiquity.csv", index=False)
    genus.to_csv(out/"accepted_genera_GloBI_moth_source.csv", index=False)
    bins.to_csv(out/"prevalence_bins_moth_documentation_rate.csv", index=False)
    report = {
      "schema": SCHEMA, "source_original_photo_IDs": 93209351,
      "matched_source_islands": 562, "matched_site_years": len(sites),
      "accepted_species_site_years": len(data),
      "unique_accepted_species": len(species_table),
      "recorded_genera": len(genus), "mixed_genera": len(mixed),
      "top_one_percent_species": len(top1),
      "top_one_percent_species_site_share": float(data.species.isin(top1).mean()),
      "top_one_percent_moth_source_species": int(
           species_table.head(top1_n).moth.sum()),
      "models": models, "exact_common_species_contributions": contribution,
      "records_are_not_real_nocturnal_pollination_or_native_flora": True,
      "observational_photo_and_GloBI_source_bias_not_identified_causally": True,
      "models_posthoc_exploratory_no_FDR_claim": True
    }
    (out/"real_93m_genus_and_prevalence_bias_result.json").write_text(
        json.dumps(report, indent=2, sort_keys=True)+"\n")
    print("ALAN_REAL_GENUS_AND_PREVALENCE_BIAS_RESULT",
          json.dumps({"model_moth": models["moth_genus_FE"],
                      "model_drop1": models["moth_genus_FE_drop_top1pct_species"],
                      "commonness_contributions": contribution},
                     sort_keys=True), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for opt in ("sites", "globi", "groups", "out"):
        p.add_argument("--"+opt, type=Path, required=True)
    a = p.parse_args()
    run(a.sites, a.globi, a.groups, a.out)


if __name__ == "__main__":
    main()
