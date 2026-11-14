#!/usr/bin/env python3
"""Actual corrected Dryad v4 frame-level moth activity experiment source audit.

DOI 10.5061/dryad.jq2bvq8ps, Proc R Soc B (2026) 10.1098/rspb.2025.2704.
Download only official Dryad /api/v2/files download links and verify source
sha256 and exact size. Count unique animals vs frames and experiment-treatment
support; inspect original R movement rules. This stage does *not* compute an
activity effect or call a source picture a pollen visit, field population or
sample of independent video frames.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,urllib.request,urllib.error
from pathlib import Path
from collections import Counter, defaultdict

FILES={
 "Activity_MainData.csv":("4507228",77107372,"e08798c04d5708e5ea51e0ceda8f4ea64f192b38656726c1f7b9a94305f18c04"),
 "Analysis_Script.R":("4507224",79730,"f33cd8af210527dbd98d8c278731f53fea0fefc01485ad0945f925478274586e"),
 "Species_size.csv":("4507222",764,"c7b872843166acbf55ca04df3662592f9bac7f5127250ea80730539b6d761e6e")
}
BASE="https://datadryad.org/api/v2/files/{}/download"
def require(a,msg):
    if not a:raise ValueError("ALAN_DRYAD23_SOURCE: "+msg)
def get(name):
    token,n,sha=FILES[name]
    request=urllib.request.Request(BASE.format(token),headers={
        "User-Agent":"ALAN-original-Dryad-moth-spectrum-source/1.0"})
    with urllib.request.urlopen(request,timeout=90) as f:
        data=f.read(n+5)
    require(len(data)==n,"Official Dryad corrected source size changed "+name)
    require(hashlib.sha256(data).hexdigest()==sha,
            "Official Dryad corrected source SHA256 changed "+name)
    return data
def median(x):
    z=sorted(x);n=len(z)
    return float(z[n//2]) if n%2 else float((z[n//2-1]+z[n//2])/2)
def audit(out):
    out=Path(out)
    require(not out.exists(),"Never overwrite existing original audit")
    try:
        original_code_bytes=get("Analysis_Script.R")
    except urllib.error.HTTPError as e:
        if e.code != 401:raise
        result={
          "schema":"ALAN_Dryad2025_v4_source_main_23moth_motion_frames_original_audit_v1",
          "status":"ORIGINAL_FILE_METADATA_VERIFIED_BYTES_BLOCKED_DRYAD_AUTH",
          "dataset_DOI":"10.5061/dryad.jq2bvq8ps",
          "source_publication_DOI":"10.1098/rspb.2025.2704",
          "source_version":4,
          "original_data_sha256_published":FILES["Activity_MainData.csv"][2],
          "original_R_code_sha256_published":FILES["Analysis_Script.R"][2],
          "Dryad_public_download_endpoint":"https://datadryad.org/api/v2/files/4507224/download",
          "HTTP_status":e.code,
          "requires_authentication_without_assuming_data_absent":True,
          "original_data_rows_recomputed":False,
          "published_paper_reported_moths":843,
          "published_paper_reported_species":23,
          "published_paper_reported_families":3,
          "published_paper_reported_10lux_white_LED_activity_suppression_pct":85,
          "original_image_frames_not_treated_independent":True,
          "not_claimed_new_species_sensitivity_coefficient":True
        }
        out.parent.mkdir(parents=True,exist_ok=True)
        out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
        print("ALAN_REAL_DRYAD_23SP_AUTH_GATE",json.dumps(result),flush=True)
        return result
    code=original_code_bytes.decode("utf8-sig",errors="replace")
    origin_code_lines=code.splitlines()
    hints=[]
    for i,x in enumerate(origin_code_lines):
        if re.search(r"mov(e|ing|ement)|activity|distance|shift",x,re.I):
            if not x.strip().startswith("#") or "movement" in x.lower():
                hints.append({"line":i+1,"source_original_R_code_line_excerpt":x[:220]})
    # Store excerpts of source R analysis only, not reupload full original
    # copyrighted source; direct source DOI is supplied in receipt.
    hints=hints[:95]
    size=get("Species_size.csv").decode("utf8-sig",errors="replace")
    taxa=list(csv.DictReader(io.StringIO(size)))
    raw=get("Activity_MainData.csv").decode("utf8-sig",errors="replace")
    source=csv.DictReader(io.StringIO(raw))
    columns=list(source.fieldnames or [])
    required={"Moth_ID","Species","Family","Light_type","Light_intensity",
              "Treatment","CaptureMethod","Time","Xcm","Ycm"}
    require(required <= set(columns),
            "Original Dryad activity source columns missing: "+str(sorted(required-set(columns))))
    per=defaultdict(lambda:{"rows":0,"species":set(),"family":set(),
                            "treatment":set(),"intensity":set(),"lighttype":set(),
                            "capture":set(),"date":set()})
    groups=Counter()
    total=0;bad=0
    for r in source:
        total+=1
        moth=str(r.get("Moth_ID") or "").strip()
        if not moth or moth=="NA":
            bad+=1
            continue
        z=per[moth]
        z["rows"]+=1
        for k,src in (("species","Species"),("family","Family"),
                      ("treatment","Treatment"),("intensity","Light_intensity"),
                      ("lighttype","Light_type"),("capture","CaptureMethod"),
                      ("date","Date")):
            z[k].add(str(r.get(src) or "").strip())
        groups[(str(r.get("Light_type") or ""),
                str(r.get("Light_intensity") or ""))]+=1
    single=lambda k:sum(len(x[k])==1 for x in per.values())
    multi={k:sum(len(x[k])>1 for x in per.values()) for k in (
             "species","family","treatment","intensity","lighttype","capture","date")}
    counts={"frames":total,"unique_nonblank_moth_IDs":len(per),
      "source_rows_missing_moth_id":bad,
      "distinct_recorded_plant_pollination_events":0,
      "distinct_moth_species":len(set.union(*(z["species"] for z in per.values()))),
      "distinct_moth_families":sorted(set.union(*(z["family"] for z in per.values()))),
      "distinct_light_types":sorted(set.union(*(z["lighttype"] for z in per.values()))),
      "distinct_illuminance_levels_original":sorted(set.union(*(z["intensity"] for z in per.values()))),
      "distinct_treatments":sorted(set.union(*(z["treatment"] for z in per.values()))),
      "capture_method_labels":sorted(set.union(*(z["capture"] for z in per.values()))),
      "median_frame_rows_per_mothID":median([z["rows"] for z in per.values()]),
      "range_frames_per_mothID":[min(z["rows"] for z in per.values()),
                                 max(z["rows"] for z in per.values())],
      "cross_design_mothIDs_with_multiple_source_levels":multi,
      "source_mothIDs_with_single_species_label":single("species"),
      "species_morphology_file_taxa":len(taxa),
      "frame_rows_by_source_lighttype_lux":[{"Light_type":k[0],
                "Light_intensity":k[1],"frame_rows":int(v)}
               for k,v in sorted(groups.items())],
      "source_columns":columns}
    require(total>100_000 and len(per)>=700 and len(per)<=1100,
            "Original 843-moth frame source grain unexpectedly changed")
    require(counts["distinct_moth_species"]==23,
            "Original published 23 species not present in actual source")
    require("Erebidae" in counts["distinct_moth_families"] and
            "Noctuidae" in counts["distinct_moth_families"],
            "Moth source families not verified")
    result={"schema":"ALAN_Dryad2025_v4_source_main_23moth_motion_frames_original_audit_v1",
       "status":"CORRECTED_ORIGINAL_DRYAD_FRAMES_INDIVIDUAL_GRAIN_VERIFIED_NO_ACTIVITY_EFFECT_ESTIMATED",
       "source":{"dataset_DOI":"10.5061/dryad.jq2bvq8ps",
         "primary_paper_DOI":"10.1098/rspb.2025.2704",
         "Dryad_corrected_public_version":4,
         "corrected_original_source_SHA256_main":FILES["Activity_MainData.csv"][2],
         "corrected_original_source_SHA256_R_script":FILES["Analysis_Script.R"][2],
         "original_source_date":"2025-12-11"},
       "source_genuine_data_counts":counts,
       "published_R_code_relevant_movement_excerpt":hints,
       "guardrails":{
         "frame_rows_not_independent_individuals":True,
         "same_moth_across_multiple_treatments_must_not_be_assumed_independent":True,
         "original_R_code_movement_definition_must_be_replicated_before_activity_effect":True,
         "natural_light_radiance_vs_lux_not_satellite_PCNL":True,
         "no_direct_flower_visit_or_pollen_deposition_measured":True,
         "no_island_mainland_predator_or_phenology_in_this_experiment":True,
         "all_original_frame_rows_not_exported":True
       }}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_REAL_DRYAD_23SP_MAIN_SOURCE",json.dumps(result,ensure_ascii=False),flush=True)
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();audit(a.out)
