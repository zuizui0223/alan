#!/usr/bin/env python3
"""Regional plant-range sensitivity for ALAN, using original WCVP TDWG3 evidence.

This joins the SAME 93,209,351-photo signed WCVP/GloBI/GSHHG island flora
source to Kew WCVP (official botanical-country / TDWG3 native-introduced
status), with each source observation year's prior 1992-2024 PCNL.

Important biological boundary: 'native_at_TDWG3' means a botanical-country
native range annotation. It does NOT verify an individual photographed plant
is wild, or is indigenous to a specific small offshore island.
Neither 'unknown' nor 'native' may be silently equated to wild-native.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import make_valid
from shapely.geometry import GeometryCollection, Polygon, MultiPolygon

from scripts.analysis import alan_real_genus_prevalence_bias_20261110 as prior
from scripts.analysis import alan_full93m_moth_PCNLPrior_two_way_FE_20261109 as model

WCVP_DIST_SHA = "9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3"
WCVP_SOURCE_SHA = "dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2"
TDWG_COMMIT = "52da7828aba9d461dd133c27b3bd7a4407161f54"
TDWG_GIT_BLOB = "91104e5159e31f88154833a51d3b0d1c9271083f"
KEY = prior.KEY
LIGHT = prior.LIGHT
STATUS = ("native_at_TDWG3", "introduced_at_TDWG3")
SCHEMA = "ALAN_genuine93m_WCVP_TDWG3_REGIONAL_status_PCNL_source_v1"


def require(ok, msg):
    if not ok:
        raise ValueError("ALAN_TDWG3_WCVP_REGIONAL: " + str(msg))


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(4 << 20), b""):
            h.update(part)
    return h.hexdigest()


def _polygon_piece(geom):
    if isinstance(geom, (Polygon, MultiPolygon)):
        return geom
    if isinstance(geom, GeometryCollection):
        from shapely import union_all
        parts = [_polygon_piece(z) for z in geom.geoms]
        parts = [z for z in parts if z is not None and not z.is_empty]
        return union_all(parts) if parts else None
    return None


def polygons(path):
    p = gpd.read_file(path)
    options = [k for k in p.columns if k.casefold() == "level3_cod"]
    require(len(options) == 1 and p.crs is not None and p.crs.to_epsg() == 4326,
            "Not original EPSG:4326 TDWG level3 polygons")
    p = p[[options[0], "geometry"]].rename(columns={options[0]: "tdwg"})
    p["tdwg"] = p.tdwg.astype(str).str.strip()
    require(p.tdwg.str.fullmatch("[A-Z]{3}").all() and len(p)>=300,
            "Original TDWG3 code universe incomplete")
    require(p.geometry.notna().all() and (~p.geometry.is_empty).all(),
            "TDWG polygons missing or empty")
    needs = ~p.geometry.is_valid
    for ix in p.index[needs]:
        polygon = _polygon_piece(make_valid(p.at[ix, "geometry"]))
        require(polygon is not None and not polygon.is_empty and polygon.is_valid,
                "Unable to repair original TDWG polygon ring")
        p.at[ix, "geometry"] = polygon
    return p


def overlay_sites(sites, boundaries):
    """One real source 0.05-degree photo site can belong to zero or one TDWG3."""
    q = sites[KEY + ["lat", "lon"]].copy().reset_index(drop=True)
    points = gpd.GeoDataFrame(
        {"row_id": np.arange(len(q),dtype=int)},
        geometry=gpd.points_from_xy(q.lon, q.lat), crs="EPSG:4326")
    hits = gpd.sjoin(points, boundaries[["tdwg", "geometry"]],
                     how="left", predicate="within")
    found = hits[["row_id", "tdwg"]].dropna().drop_duplicates()
    n = found.groupby("row_id",observed=True).tdwg.agg(["nunique","first"])
    q["tdwg_level3_code"] = ""
    q["tdwg_region_state"] = "outside_or_boundary"
    unique = n.loc[n["nunique"].eq(1)]
    q.loc[unique.index, "tdwg_level3_code"] = unique["first"].to_numpy()
    q.loc[unique.index, "tdwg_region_state"] = "unique_pinned_TDWG3"
    q.loc[n.index[n["nunique"].gt(1)], "tdwg_region_state"] = "overlapping_TDWG3_unknown"
    return q[KEY+["tdwg_level3_code","tdwg_region_state"]]


def regional_table(path):
    require(sha256(path)==WCVP_DIST_SHA, "Official WCVP distribution SHA changed")
    names = ["wcvp_accepted_plant_name_id", "wcvp_accepted_name",
             "tdwg_level3_code", "native_introduced_state"]
    d = pd.read_csv(path,compression="gzip",dtype=str,keep_default_na=False,usecols=names)
    require(len(d)==657091 and d.wcvp_accepted_plant_name_id.nunique()==73232,
            "Known 657,091 WCVP source country×taxon records drifted")
    require(not d.duplicated(["wcvp_accepted_name","tdwg_level3_code"]).any(),
            "One species and region with contradictory WCVP country-native state")
    require(set(d.native_introduced_state).issubset(
            {"native_at_TDWG3","introduced_at_TDWG3","extinct","doubtful",
             "unknown_flag","contradictory_TDWG3_records"}),
            "Unrecognized original WCVP range category")
    return d.rename(columns={"wcvp_accepted_name":"species"})


def source_species(groups, sites):
    """True WCVP accepted species×site-year, photo IDs deduped across 16 groups."""
    frames=[]
    cols = ["island","latbin","lonbin","year","species",
            "moth_positive_photo_ids","butterfly_positive_photo_ids",
            "moth_pollinates_photo_ids","night_pollen_photo_ids"]
    for i in range(16):
        file=Path(groups)/f"allflora-island-group-{i:02d}.sqlite"
        receipt=file.with_suffix(".json")
        require(file.is_file() and receipt.is_file(),
                f"Missing verified 93m original source group {i}")
        audit=json.loads(receipt.read_text())
        require(audit.get("group")==i and
                audit.get("schema")=="ALAN_full93m_genuine_WCVP_allflora_moth_evidence_island_source_group_v1"
                and audit.get("counts",{}).get("completed_verified_original_tiers")==14,
                "Original 14-tier photo source signature missing")
        with sqlite3.connect(file) as con:
            d=pd.read_sql_query("SELECT "+",".join(cols)+" FROM flora",con)
        d["island"]=d.island.astype(str)
        d = d.merge(sites[KEY],on=KEY,how="inner",validate="many_to_one")
        frames.append(d)
    x=pd.concat(frames,ignore_index=True)
    f=("moth_positive_photo_ids","butterfly_positive_photo_ids",
       "moth_pollinates_photo_ids","night_pollen_photo_ids")
    for k in f:
        x[k]=x[k].gt(0)
    duplicates=x.duplicated(KEY+["species"],keep=False)
    if duplicates.any():
        agree=x.loc[duplicates].groupby(KEY+["species"],observed=True)[list(f)].nunique()
        require((agree.to_numpy()<=1).all(),
                "Original groups disagree on accepted species positive guild evidence")
    x=x.drop_duplicates(KEY+["species"])
    counts=x.groupby(KEY,observed=True).agg(
        n_species=("species","size"),
        moth_species=("moth_positive_photo_ids","sum")).reset_index()
    checked=sites[KEY+["all_accepted_recorded_species",
                       "moth_visitor_positive_recorded_species"]].merge(
        counts,on=KEY,how="left",validate="one_to_one")
    require(len(x)==2012925 and
            checked.n_species.eq(checked.all_accepted_recorded_species).all() and
            checked.moth_species.eq(checked.moth_visitor_positive_recorded_species).all(),
            "93m WCVP taxon×site-year source does not conserve original data")
    return x


def summarise_status(species, places, regional, sites, *, return_species=False):
    x = species.merge(places[KEY+["tdwg_level3_code","tdwg_region_state"]],
                      on=KEY,how="left",validate="many_to_one")
    d=regional[["species","tdwg_level3_code","native_introduced_state"]]
    x=x.merge(d,on=["species","tdwg_level3_code"],how="left",
              validate="many_to_one")
    x["range_class"]=x.native_introduced_state.fillna("unknown")
    x.loc[~x.tdwg_region_state.eq("unique_pinned_TDWG3"),"range_class"]="unknown"
    x.loc[~x.range_class.isin(STATUS),"range_class"]="unknown"
    # This is the region-wide STATUS of the SPECIES, NOT a photo's wild status.
    x["moth"]=x.moth_positive_photo_ids.astype(int)
    x["butterfly"]=x.butterfly_positive_photo_ids.astype(int)
    x["night_pollen"]=x.night_pollen_photo_ids.astype(int)
    grp=x.groupby(KEY+["range_class"],observed=True).agg(
        recorded_species=("species","size"),
        recorded_moth=("moth","sum"),
        recorded_butterfly=("butterfly","sum"),
        recorded_night_pollen=("night_pollen","sum")
    ).reset_index()
    p=grp.pivot(index=KEY,columns="range_class",
                values=["recorded_species","recorded_moth",
                        "recorded_butterfly","recorded_night_pollen"])
    p.columns=["__".join(k) for k in p.columns]
    p=p.reset_index()
    r=sites.merge(p,on=KEY,how="left",validate="one_to_one")
    cats=list(STATUS)+["unknown"]
    for status in cats:
        for count in ("recorded_species","recorded_moth",
                      "recorded_butterfly","recorded_night_pollen"):
            col=count+"__"+status
            if col not in r: r[col]=0
            r[col]=r[col].fillna(0).astype(int)
    require(r[[f"recorded_species__{j}" for j in cats]].sum(axis=1).eq(
            r.all_accepted_recorded_species).all(),
            "WCVP native/introduced/unknown species do not conserve source denominator")
    require(r[[f"recorded_moth__{j}" for j in cats]].sum(axis=1).eq(
            r.moth_visitor_positive_recorded_species).all(),
            "GloBI moth records not conserved in WCVP native/introduced/unknown states")
    return (r,x) if return_species else r


def model_group(panel, status, *, min_species=3, min_islands=30):
    den="recorded_species__"+status
    num="recorded_moth__"+status
    x=panel.loc[panel[den].ge(min_species)].copy()
    # A region-classified site with one sampled site cannot identify a within-
    # island light contrast. Require same-year spatially distinct light values.
    mix=x.groupby(["island","year"],observed=True)[LIGHT].agg(["size","min","max"])
    valid=mix.loc[mix["size"].ge(2)&(mix["max"]-mix["min"]).ge(.02)]
    keys=set(valid.reset_index().island)
    x=x.loc[x.island.isin(keys)].copy()
    if len(x)<100 or x.island.nunique()<min_islands:
        return {"status":"not_estimable","eligible_sites":len(x),
                "eligible_islands":int(x.island.nunique()),
                "reason":"too few spatially supported islands after WCVP region-state filter"}
    x["frac_pp"]=100*x[num]/x[den]
    x["log_obs"]=np.log1p(x.all_accepted_observation_uuids)
    x["log_species"]=np.log1p(x[den])
    # Existing exploratory model: island FE + photo-year FE, island equal.
    base_fit=model.fit(x,"island_equal")
    site_fit=model.fit(x,"site_equal")
    species_fit=model.fit(x,"island_equal",True)
    # More conservative contemporaneous island×year FE: removes all island-
    # specific annual temporal shifts, including sampling pulses.
    iy=x.island.astype(str).str.cat(x.year.astype(str),sep="|")
    count=x.groupby("island",observed=True).size()
    weights=1/x.island.map(count).to_numpy(float)
    # Slope calculated in PP, rather than log odds; avoid unobserved zeros.
    from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import fixed_effects_regression
    fit_iy=fixed_effects_regression(
        x[num].to_numpy(float)/x[den].to_numpy(float),
        x[LIGHT].to_numpy(float),x.log_obs.to_numpy(float),
        weights,x.island.to_numpy(),[iy.to_numpy()])
    for label,f in (("island_plus_year_FE",base_fit),
                    ("site_equal_island_plus_year_FE",site_fit),
                    ("island_plus_year_FE_richness_adjusted",species_fit)):
        require(f.get("status")=="estimated_exploratory_observational",
                "Original WCVP regional %s model rank failed" % label)
    return {
      "status":"estimated_exploratory_observational",
      "eligible_sites":len(x),"eligible_islands":int(x.island.nunique()),
      "accepted_species_site_records":int(x[den].sum()),
      "moth_positive_species_site_records":int(x[num].sum()),
      "island_equal_island_plus_year_FE":base_fit["terms"][LIGHT],
      "site_equal_island_plus_year_FE":site_fit["terms"][LIGHT],
      "island_equal_plus_species":species_fit["terms"][LIGHT],
      "island_equal_island_x_year_FE":fit_iy
    }


def native_genus_bias(source_species_status, original_sites):
    """Same WCVP-region-native source subset, species-level genus/ubiquity guards.

    The top 1% species are ranked only by observed source geographic frequency,
    never by ALAN coefficient or flower colour. A regional-native class is
    not proof that an iNaturalist photographed individual is wild.
    """
    from scripts.analysis import alan_real_genus_prevalence_bias_20261110 as genus
    d=source_species_status.copy()
    require(len(d)==2012925, "WCVP original source accepted species-site population changed")
    d["genus"]=d.species.astype(str).str.extract(
        r"^([A-Z][a-zA-Z-]+)\\b",expand=False).fillna("UNKNOWN_GENUS")
    d["moth_minus_butterfly"]=d.moth.astype(int)-d.butterfly.astype(int)
    prevalence=d.groupby("species",observed=True).size().sort_values(ascending=False)
    require(len(prevalence)==14330, "Original 14,330 accepted WCVP source species changed")
    topn=int(np.ceil(len(prevalence)*.01))
    require(topn==144, "Prespecified 1% WCVP species source size changed")
    top_species=set(prevalence.head(topn).index)
    # No unsupported geographical absence: all species are actually photographed.
    regional=d.loc[d.range_class.eq("native_at_TDWG3")].copy()

    def supported(table):
        count=table.groupby(KEY,observed=True).size().reset_index(name="n_recorded")
        src=original_sites.merge(count,on=KEY,how="inner",validate="one_to_one")
        src=src.loc[src.n_recorded.ge(3)].copy()
        if src.empty:return src
        bounds=src.groupby(["island","year"],observed=True)[LIGHT].agg(
            ["size","min","max"])
        ok=set(bounds.loc[bounds["size"].ge(2)&
                (bounds["max"]-bounds["min"]).ge(.02)].reset_index().island)
        return src.loc[src.island.isin(ok)].copy()

    def fitted(table,support,response,genus_fe):
        if len(support)<100 or support.island.nunique()<30:
            return {"status":"not_estimable","n_sites":len(support),
                    "n_islands":int(support.island.nunique())}
        return genus.species_model(table,support,response=response,
                                   genus_fe=genus_fe,weight="island_equal")

    support=supported(regional)
    baseline=fitted(regional,support,"moth",False)
    require(abs(baseline["beta_pp"]-0.9155117901610345)<1e-6,
            "WCVP-regional native source model must recover original island×year FE")
    models={
       "native_moth_without_genus_FE":baseline,
       "native_moth_with_genus_FE":fitted(regional,support,"moth",True),
       "native_butterfly_with_genus_FE":fitted(regional,support,"butterfly",True),
       "native_moth_minus_butterfly_genus_FE_direct":fitted(
           regional,support,"moth_minus_butterfly",True)
    }
    if (models["native_moth_with_genus_FE"].get("beta_pp") is not None
        and models["native_butterfly_with_genus_FE"].get("beta_pp") is not None):
        require(abs(models["native_moth_with_genus_FE"]["beta_pp"]-
                    models["native_butterfly_with_genus_FE"]["beta_pp"]-
                    models["native_moth_minus_butterfly_genus_FE_direct"]["beta_pp"])<1e-7,
                "Genus adjusted same source guild difference did not close")
    subset=regional.loc[~regional.species.isin(top_species)].copy()
    s2=supported(subset)
    models["native_genus_FE_excluding_144_most_observed_global_source_species"]=(
        fitted(subset,s2,"moth",True))
    models["native_butterfly_genus_FE_excluding_144_most_observed_global_species"]=(
        fitted(subset,s2,"butterfly",True))
    return {
      "original_source_species":int(len(prevalence)),
      "global_top_one_percent_source_species":topn,
      "global_top_one_percent_region_native_species_site_fraction":float(
          regional.species.isin(top_species).mean()),
      "region_native_source_species_site_rows":len(regional),
      "region_native_sites_with_min3_species_and_spatial_light":len(support),
      "region_native_islands_with_spatial_light":int(support.island.nunique()),
      "models":models,
      "source_common_species_selected_before_ALAN_outcome":True,
      "native_WCVP_region_NOT_wild_island_native":True,
      "all_models_exploratory_after_bias_discovery":True
    }


def run(sites, groups, wcvp, geojson, out):
    s=prior.load_sites(sites)
    assert len(s)==77647
    native=regional_table(wcvp)
    border=polygons(geojson)
    sites_geo=overlay_sites(s,border)
    species=source_species(groups,s)
    joined,species_status=summarise_status(species,sites_geo,native,s,return_species=True)
    native_genus_audit=native_genus_bias(species_status,s)
    cats=list(STATUS)+["unknown"]
    counts={k:{
       "recorded_species_site_years":int(joined[f"recorded_species__{k}"].sum()),
       "GloBI_moth_positive_species_site_years":int(joined[f"recorded_moth__{k}"].sum()),
       "n_sites_with_any_recorded_species":int(joined[f"recorded_species__{k}"].gt(0).sum())
    } for k in cats}
    total=joined.all_accepted_recorded_species.sum()
    require(sum(v["recorded_species_site_years"] for v in counts.values())==total,
            "WCVP regional source species totals not conserved")
    versions={
       "WCVP_regional_native_proxy":model_group(joined,"native_at_TDWG3"),
       "WCVP_regional_introduced_proxy":model_group(joined,"introduced_at_TDWG3"),
       "WCVP_unknown_or_unmatched":model_group(joined,"unknown")
    }
    # Same geographic sample pair: require >=3 documented native and
    # introduced WCVP-status species at a source cell/year; missingness isn't 0.
    paired=joined.loc[joined["recorded_species__native_at_TDWG3"].ge(3)&
                      joined["recorded_species__introduced_at_TDWG3"].ge(3)]
    if len(paired)>=100 and paired.island.nunique()>=30:
        paired=paired.copy()
        paired["frac_pp"]=100*(paired["recorded_moth__native_at_TDWG3"]/
                               paired["recorded_species__native_at_TDWG3"]-
                               paired["recorded_moth__introduced_at_TDWG3"]/
                               paired["recorded_species__introduced_at_TDWG3"])
        paired["log_obs"]=np.log1p(paired.all_accepted_observation_uuids)
        paired["log_species"]=np.log1p(paired.all_accepted_recorded_species)
        # Keep only same-year/within-island contrast support.
        mix=paired.groupby(["island","year"],observed=True)[LIGHT].agg(["size","min","max"])
        keep=set(mix.loc[mix["size"].ge(2)&
                (mix["max"]-mix["min"]).ge(.02)].reset_index().island)
        paired=paired.loc[paired.island.isin(keep)].copy()
        versions["native_minus_introduced_same_sites"]={
          "eligible_site_years":len(paired),"eligible_islands":int(paired.island.nunique())}
        if len(paired)>=100 and paired.island.nunique()>=30:
            x=model.fit(paired,"island_equal")
            versions["native_minus_introduced_same_sites"]["result"]=(
                x["terms"][LIGHT] if x.get("status")==
                "estimated_exploratory_observational" else x)
        else:
            versions["native_minus_introduced_same_sites"]["result"]="not_estimable"
    else:
        versions["native_minus_introduced_same_sites"]={
            "eligible_site_years":len(paired),
            "eligible_islands":int(paired.island.nunique()),
            "result":"not_estimable"}
    report={
      "schema":SCHEMA,
      "original_93m_source_observation_scope":93209351,
      "source_2019_2025_PCNL_site_years":len(s),
      "source_islands":int(s.island.nunique()),
      "original_site_accepted_species_total":int(total),
      "unique_TDWG3_assigned_sites":int(sites_geo.tdwg_region_state.eq("unique_pinned_TDWG3").sum()),
      "WCVP_status_summary":counts,
      "models":versions,
      "native_WCVP_genus_ubiquitous_taxa_source_audit":native_genus_audit,
      "WCVP_status_is_botanical_country_not_island_native":True,
      "photographed_individual_wild_or_garden_status_unknown":True,
      "missing_distribution_unknown_not_introduced":True,
      "observed_source_GloBI_moth_record_not_effective_moth_pollination":True,
      "new_estimates_exploratory_not_causal":True,
      "TDWG_source_commit":TDWG_COMMIT,
      "TDWG_source_git_blob":TDWG_GIT_BLOB,
      "TDWG_source_SHA256":sha256(geojson),
      "WCVP_distribution_SHA256":WCVP_DIST_SHA
    }
    folder=Path(out)
    require(not folder.exists(),"Refuse to replace previous real WCVP overlay")
    folder.mkdir(parents=True)
    joined.to_csv(folder/"ALAN_real_PCNL_WCVP_TDWG3_regional_status_island_site_year.csv.gz",
                  index=False,compression="gzip")
    (folder/"ALAN_real_PCNL_WCVP_TDWG3_origin_sensitivity_20261110.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_REAL_93M_WCVP_TDWG3_REGION_NOT_ISLAND_NATIVE_RESULT",
          json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("sites","groups","wcvp","geojson","out"):
        p.add_argument("--"+k,type=Path,required=True)
    args=p.parse_args()
    run(args.sites,args.groups,args.wcvp,args.geojson,args.out)

if __name__=="__main__":
    main()
