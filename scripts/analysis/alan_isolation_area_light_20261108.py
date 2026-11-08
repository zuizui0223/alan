#!/usr/bin/env python3
"""ALAN: nearest continental shoreline versus island area as nightlight-slope moderators.

This is an exploratory, source-verified OBSERVATIONAL comparison, not a
pollination/selection experiment. Mainland isolation is estimated as the minimum
WGS84 GEODESIC distance between GSHHG high-resolution ORIGINAL OUTER COASTLINE
VERTICES for each island and any mainland landmass defined in the original
study (>7 million km²). It is a VERTEX-BASED approximation to the continuous
coast-to-coast distance, NOT a centroid/island-capital or ferry route distance.

2019/21/23/24/25 original-photo true-ROI nine-palette, exact WCVP taxon,
same-year VIIRS & WorldPop are verified using the ALAN-native inherited loaders.
Unknown colours, photo-ineligible sites and missing years are never imputed.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import Geod, Transformer
from scipy.spatial import cKDTree
import statsmodels.formula.api as smf

from scripts.analysis import alan_island_area_light_20261108 as area
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import prepare_land
from scripts.analysis.mois_global2025_original_source_20261107 import COAST_SHA

YEARS=(2019,2021,2023,2024,2025)
EXPECTED={"strict_5p3r3s":7982,"measured_3p2r2s":20543}
BASE_AREA_BETAS={
 "strict_5p3r3s":{"cell":2.208728309637273,"island":5.151623122161728},
 "measured_3p2r2s":{"cell":0.9983139438500241,"island":3.3846627689849633}
}
BASE_COVARS="log_pop + bio5 + log_rain + elevation + log_attempts"
TERMS=("log_light","log_light:log10_area_c","log_light:log10_distance_c")
TRANSFORMER=Transformer.from_crs("EPSG:4326","EPSG:4978",always_xy=True)
GEOD=Geod(ellps="WGS84")

def load_signed_land(coast):
    coast=Path(coast)
    receipt=json.loads((coast/"source_receipt.json").read_text())
    panel.require(receipt.get("original_archive_sha256")==COAST_SHA,
                  "Original GSHHG 2.3.7 archive SHA is not frozen")
    for rel,digest in receipt["extracted_sha256"].items():
        panel.require(panel.sha(coast/rel)==digest,
                      "Extracted GSHHG bytes drift: "+rel)
    land=prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp"))
    land.land_id=land.land_id.astype(str)
    panel.require(land.land_id.is_unique and land.area_km2.gt(0).all(),
                  "Invalid source GSHHG landmass geometry or area")
    return land

def outer_vertices(geometry):
    """Only the outer coastline of each polygon part; omit inland holes."""
    parts=shapely.get_parts(geometry)
    rings=shapely.get_exterior_ring(parts)
    pts=shapely.get_coordinates(rings)
    panel.require(len(pts)>=4 and np.isfinite(pts).all(),
                  "Geometry has too few valid source outer-coast vertices")
    return pts

def geocentric(pts):
    xx,yy,zz=TRANSFORMER.transform(
        np.asarray(pts[:,0],dtype=float),
        np.asarray(pts[:,1],dtype=float),
        np.zeros(len(pts),dtype=float))
    coords=np.column_stack((xx,yy,zz))
    panel.require(np.isfinite(coords).all(),
                  "Nonfinite source coastline geocentric coordinates")
    return coords

def isolated_mainland_distance_km(land,source_island_ids):
    """Approximate nearest mainland coast using original WGS84 coast VERTICES.

    ECEF nearest-neighbour search is globally safe across the antimeridian.
    For each island vertex, its nearest three mainland coastal vertices are
    compared using WGS84 geodesic distance, in bounded 15k-vertex batches.
    Thus all observed source vertices are examined with no centroid proxy.
    Distances are an upper bound approximation to exact segment distance.
    """
    src_ids=sorted(set(map(str,source_island_ids)))
    classified=land.set_index("land_id")
    panel.require(len(src_ids)>0 and set(src_ids).issubset(set(classified.index)),
                  "Requested GSHHG islands absent from original polygons")
    panel.require(classified.loc[src_ids,"land_class"].eq("island").all(),
                  "A source mainland was used as an island")
    main=land.loc[land.land_class.eq("mainland")].copy()
    panel.require(len(main)>=1,"Original GSHHG mainland category is empty")
    all_main=np.concatenate([outer_vertices(v) for v in main.geometry],axis=0)
    panel.require(len(all_main)>=10,"Too few original mainland coastline vertices")
    tree=cKDTree(geocentric(all_main))
    neighbour_k=min(3,len(all_main))
    result={}
    for ident in src_ids:
        vertices=outer_vertices(classified.loc[ident,"geometry"])
        best=np.inf
        for start in range(0,len(vertices),15000):
            batch=vertices[start:start+15000]
            _,ix=tree.query(geocentric(batch),k=neighbour_k,workers=2)
            ix=np.asarray(ix,dtype=int).reshape(-1,neighbour_k)
            near=all_main[ix.reshape(-1)]
            obs=np.repeat(batch,neighbour_k,axis=0)
            _,_,dist_m=GEOD.inv(obs[:,0],obs[:,1],near[:,0],near[:,1])
            best=min(best,float(np.min(dist_m)))
        panel.require(np.isfinite(best) and best>=0,
                      "Nonfinite mainland coast-vertex distance")
        result[ident]=best/1000
    return pd.DataFrame({"land_id":list(result),
                         "distance_mainland_vertex_km":list(result.values())})

def complete_island_frame(within,distances):
    d=within.merge(distances,on="land_id",how="left",validate="many_to_one")
    panel.require(len(d)==len(within) and d.distance_mainland_vertex_km.notna().all()
                  and d.distance_mainland_vertex_km.ge(0).all(),
                  "GSHHG original island-distance join was incomplete")
    per=d.groupby("land_id",observed=True).distance_mainland_vertex_km.nunique()
    panel.require(per.eq(1).all(),"Island distance varies within fixed landmass")
    return d

def fit_light_moderation(frame,area_centre,distance_centre,model,equal_island):
    d=frame.copy()
    d["log10_area_c"]=np.log10(d.area_km2)-area_centre
    d["log10_distance_c"]=np.log10(1+d.distance_mainland_vertex_km)-distance_centre
    panel.require(d.land_id.nunique()>=10 and d.log10_area_c.std()>0 and
                  d.log10_distance_c.std()>0,
                  "Insufficient distinct islands or modifier variation")
    modifiers={
       "area_only":"log_light:log10_area_c",
       "distance_only":"log_light:log10_distance_c",
       "joint":"log_light:log10_area_c + log_light:log10_distance_c"
    }
    panel.require(model in modifiers,"Unregistered area/isolation comparison")
    formula=("white_pp ~ log_light + "+modifiers[model]+" + "+BASE_COVARS+
             " + C(land_id) + C(year)")
    if equal_island:
        weight=1/d.groupby("land_id",observed=True).land_id.transform("size").astype(float)
        panel.require(np.allclose(weight.groupby(d.land_id).sum(),1.0),
                      "Island-equal weighting not one per landmass")
        fit=smf.wls(formula,d,weights=weight).fit(
            cov_type="cluster",cov_kwds={"groups":d.land_id},use_t=True)
    else:
        fit=smf.ols(formula,d).fit(
            cov_type="cluster",cov_kwds={"groups":d.land_id},use_t=True)
    output={}
    for term in TERMS:
        if term not in fit.params:
            continue
        lo,hi=fit.conf_int().loc[term]
        vals={"beta":float(fit.params[term]),"SE":float(fit.bse[term]),
              "ci_low":float(lo),"ci_high":float(hi),"p":float(fit.pvalues[term])}
        panel.require(np.isfinite(list(vals.values())).all(),
                      "Nonfinite observed island/area/distance regression")
        output[term]=vals
    return {"terms":output,"model":formula,
            "weighting":"equal_total_per_island" if equal_island
                        else "equal_per_observed_cell_year",
            "n_islands":int(d.land_id.nunique()),"n_cell_years":len(d),
            "island_clustered_nominal_exploratory_CI":True}

def run(paths,coast,out):
    out=Path(out)
    panel.require(not out.exists(),"Refuse overwrite of original ALAN results")
    frames=tier.build_frozen_five_years(paths)
    cohorts=tier.by_year_tier(frames)
    land=load_signed_land(coast)
    meta=land[["land_id","land_class","area_km2"]]
    selected={}
    for name,expected in EXPECTED.items():
        q=cohorts[name]
        panel.require(len(q)==expected and q.target_id.is_unique and
                      set(q.year)==set(YEARS),
                      "Original real-photo/VIIRS source selection drift")
        all_island,within,groups=area.join_area_and_select(q,meta)
        selected[name]=(all_island,within,groups)
    core=selected["strict_5p3r3s"]
    panel.require(len(core[0])==710 and core[1].land_id.nunique()==44 and
                  len(core[1])==650,"Prespecified strict island source changed")
    broad=selected["measured_3p2r2s"]
    panel.require(len(broad[0])==2026 and broad[1].land_id.nunique()==95 and
                  len(broad[1])==1884,"Prespecified broad island source changed")
    ids=set(core[1].land_id)|set(broad[1].land_id)
    coast_distance=isolated_mainland_distance_km(land,ids)
    centre_area=float(np.log10(core[1].groupby("land_id").area_km2.first()).median())
    core_distance=complete_island_frame(core[1],coast_distance)
    centre_distance=float(np.log10(1+core_distance.groupby("land_id")
                                      .distance_mainland_vertex_km.first()).median())
    outcomes={}
    for name,(_,within,_) in selected.items():
        d=complete_island_frame(within,coast_distance)
        per=d.groupby("land_id").agg(area=("area_km2","first"),
                  distance=("distance_mainland_vertex_km","first"),
                  n=("target_id","size"))
        loga=np.log10(per.area)
        logd=np.log10(1+per.distance)
        rho=float(np.corrcoef(loga,logd)[0,1])
        panel.require(np.isfinite(rho),"Island area/isolation correlation unidentified")
        fits={}
        for weights,eq in (("cell_equal",False),("island_equal",True)):
            fits[weights]={}
            for kind in ("area_only","distance_only","joint"):
                fits[weights][kind]=fit_light_moderation(
                     d,centre_area,centre_distance,kind,eq)
            old=fits[weights]["area_only"]["terms"]["log_light:log10_area_c"]["beta"]
            fixed=BASE_AREA_BETAS[name]["island" if eq else "cell"]
            panel.require(abs(old-fixed)<1e-5,
                          "Original independently signed area coefficient changed")
        outcomes[name]={
            "islands":len(per),"cell_years":len(d),
            "area_isolation_log_correlation_equal_island":rho,
            "distance_mainland_vertex_km_median_equal_island":float(per.distance.median()),
            "distance_mainland_vertex_km_min":float(per.distance.min()),
            "distance_mainland_vertex_km_max":float(per.distance.max()),
            "island_area_log10_centre_strict_equal_island":centre_area,
            "isolation_log10_1plus_km_centre_strict_equal_island":centre_distance,
            "models":fits
        }
    report={
        "schema":"ALAN_real_source_GSHHG_area_vs_mainland_isolation_light_v1",
        "source_photo_years":list(YEARS),
        "distance_definition":"nearest original WGS84 outer-coastline VERTICES (not continuous polygon distance), GSHHG L1 island to GSHHG mainland >7 million km²; WGS84 geodesic nearest among 3 3D-ECEF neighbours per island vertex",
        "distance_is_approximate_vertex_upper_bound_km":True,
        "area_definition":"original GSHHG island surface area, EPSG:6933",
        "years_are_not_independent_islands":True,
        "unverified_causal_light_pollination_selection":True,
        "no_imputed_photo_colour_or_plant_absence":True,
        "new_original_images_or_ROI_inferences":0,
        "per_island_distance_values_not_publicly_exported":True,
        "results":outcomes,
        "source_SHA256":{**{k:panel.sha(v) for k,v in paths.items()},
                        "GSHHG_archive":COAST_SHA}
    }
    out.mkdir(parents=True)
    (out/"island_isolation_area_light_aggregate_receipt.json").write_text(
           json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_ACTUAL_ISOLATION_AREA_LIGHT_RESULTS",json.dumps({
        key:{
            "area_isolation_r":v["area_isolation_log_correlation_equal_island"],
            "median_isolation_km":v["distance_mainland_vertex_km_median_equal_island"],
            "n_islands":v["islands"],
            "cell_years":v["cell_years"],
            "models":{w:{
                m:{t:{"beta":r["beta"],"ci_low":r["ci_low"],
                       "ci_high":r["ci_high"],"p":r["p"]}
                   for t,r in v["models"][w][m]["terms"].items()
                   if t.startswith("log_light:")}
                for m in ("area_only","distance_only","joint")}
                for w in ("cell_equal","island_equal")}}
        for key,v in outcomes.items()},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    keys=("canonical","wcvp","geo2019","viirs2019","pop2019",
          "geo2021","viirs2021","pop2021",
          "geo2023","viirs2023","pop2023",
          "geo2024","cn2024","uc2024","uc2024receipt",
          "geo2025","geo2025receipt","viirs2025")
    for key in keys:
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--coast",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in keys},a.coast,a.out)
