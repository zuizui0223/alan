#!/usr/bin/env python3
"""Author-original Dryad 2013 wetland emergence taxon audit, non-promoting.

Check if Sullivan2019 aquatic emergent trap Chironomidae can be source-
matched to Ogbeide2026 genuinely measured Chironomidae emergence.
Only three habitat/group blocks are independent, NOT 15 stations/4 weeks.
Do not redistribute original author rows or infer missing taxa as zero.
"""
import argparse,csv,hashlib,io,json,zipfile,math
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlparse

URL="https://datadryad.org/api/v2/versions/8787/download"
FILENAME="Sullivan_EAP18-0555_WetlandFieldExperiment_2013.txt"
DIGEST="385444ab318ec4b15edd3161c740aa7f4b2d9b973fe4c2ddf6c29cc56cdec484"
ARMS=("DRK","LEDI","LEDH","HPSI","HPSH")

def get_original():
    with urlopen(Request(URL,headers={"User-Agent":"ALAN-Verified-Source-Study-Audit/2026-11-24"}),timeout=25) as f:
        final=urlparse(f.geturl())
        if final.scheme!="https" or not final.hostname or not (
            final.hostname in {"datadryad.org","www.datadryad.org"} or
            final.hostname.endswith(".amazonaws.com")
        ):
            raise ValueError("STOP source redirected outside official Dryad HTTPS/S3 delivery")
        z=f.read(1500000)
    with zipfile.ZipFile(io.BytesIO(z)) as f:
        entries=[x for x in f.infolist() if Path(x.filename).name==FILENAME]
        if len(entries)!=1 or entries[0].file_size>400000:
            raise ValueError("STOP expected unique 2013 author study table")
        raw=f.read(entries[0])
    if hashlib.sha256(raw).hexdigest()!=DIGEST:
        raise ValueError("STOP original TSV SHA256 changed")
    return raw

def analyze(raw):
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t")
    required={"Station","Treatment","Group","Sample","Date","Taxa","Count"}
    if not required.issubset(reader.fieldnames or ()):
        raise ValueError("STOP unexpected source column schema")
    src=list(reader)
    if len(src)!=1161:raise ValueError("STOP original rows changed")
    aquatic=[r for r in src if r["Sample"].strip()=="AquaticEmergent"]
    if len(aquatic)!=198:raise ValueError("STOP original aquatic source rows changed")
    tax=Counter();per_arm=defaultdict(Counter);per_block=defaultdict(Counter)
    dates=defaultdict(set);units=defaultdict(set)
    for r in aquatic:
        grp,tr=r["Group"].strip(),r["Treatment"].strip()
        if grp not in ("1","2","3") or tr not in ARMS:
            raise ValueError("STOP unexpected block or treatment")
        day,site,name=r["Date"].strip(),r["Station"].strip(),r["Taxa"].strip()
        if not day or not site or not name or not site.endswith("-"+grp):
            raise ValueError("STOP source event/taxon key broken")
        n=float(r["Count"].strip())
        if not math.isfinite(n) or n<0 or not n.is_integer():
            raise ValueError("STOP source count not nonnegative integer")
        n=int(n)
        tax[name]+=n
        per_arm[tr][name]+=n
        per_block[grp,tr][name]+=n
        units[grp,tr].add(site);dates[grp,tr].add(day)
    if sum(tax.values())!=5372 or len(units)!=15:
        raise ValueError("STOP original aquatic counts/stations changed")
    if any(len(units[g,t])!=1 or len(dates[g,t])!=4
           for g in ("1","2","3") for t in ARMS):
        raise ValueError("STOP no four repeated dates within three source blocks")
    chiro=sorted(t for t in tax if "chironom" in t.lower())
    ambiguous=sorted(t for t in tax if any(q in t.lower() for q in
        ("diptera","midge","unidentif","unknown","other")))
    return {
        "status":"SOURCE_LABEL_RETRIEVED_NO_META_EFFECT",
        "publication_doi":"10.1002/eap.1821",
        "original_dataset_doi":"10.5061/dryad.gb9534p",
        "original_2013_author_TSV_sha256":DIGEST,
        "total_original_records":1161,
        "aquatic_emergence_source_rows":198,
        "independent_habitat_blocks":3,
        "original_stations":15,
        "total_aquatic_catches":5372,
        "distinct_source_taxon_labels":len(tax),
        "source_taxon_count_totals":dict(tax.most_common()),
        "chironomid_source_labels":chiro,
        "possibly_ambiguous_dipteran_taxon_labels":ambiguous,
        "chironomid_count_by_treatment":{tr:sum(per_arm[tr][t] for t in chiro) for tr in ARMS},
        "source_block_treatment_counts":[{
            "group":g,"light":tr,
            "all_aquatic_capture":sum(per_block[g,tr].values()),
            "explicit_chironomid_recorded_capture":sum(per_block[g,tr][t] for t in chiro),
            "four_original_dates":len(dates[g,tr]),
        } for g in ("1","2","3") for tr in ARMS],
        "meta_k_added":0,
        "limit":"Taxon labels must be adjudicated. No claim of source family completeness, sample-absence zeros or matching catch sampling effort to Ogbeide; both studies differ in site and lighting context.",
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path,default=Path("results/sullivan2019_chironomid_source_gate_20261124.json"))
    p=a.parse_args()
    res=analyze(get_original())
    p.output.parent.mkdir(parents=True,exist_ok=True)
    p.output.write_text(json.dumps(res,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("SULLIVAN_TAXON_SOURCE_GATE",json.dumps({
        "source_original_verified":True,
        "source_taxa":res["source_taxon_count_totals"],
        "Chironomidae_labels":res["chironomid_source_labels"],
        "source_chiro_counts_by_treatment":res["chironomid_count_by_treatment"],
        "physical_n":3,
        "meta_k_added":0
    },ensure_ascii=False))
if __name__=="__main__":
    main()
