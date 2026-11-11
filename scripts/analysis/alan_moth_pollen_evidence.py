"""Partition the frozen 78-region response by externally measured moth pollen evidence."""
import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path

from docx import Document
import numpy as np
import pandas as pd

from alan_gecko_primary_fiveyear import ID, primary
from alan_globi_v011_same_fiveyear import taxonomy
from alan_globi_crossland_light_robust_20261111 import fit
from mois_alan_atwater_moth_plant_source_join_20261106 import TABLE2


def pollen_records(folder):
    rows = []
    def add(names, source, measure):
        rows.extend((name, source, measure) for name in names)

    walton = pd.read_excel(folder/"rsbl20190877_si_001.xlsx", sheet_name=2, header=None).iloc[4:-1]
    assert len(walton.columns) == 12
    moths = walton.iloc[:,5:12].apply(pd.to_numeric,errors="raise").sum(axis=1)
    names = walton.loc[moths.gt(0),0].str.strip().str.rstrip("*").replace({"Iris psuedacorus":"Iris pseudacorus"})
    add(names, "10.6084/m9.figshare.c.4955030", "body pollen")
    with zipfile.ZipFile(folder/"ceh.zip") as archive:
        ceh = pd.read_csv(archive.open("data/F4_pollen.csv"))
    nz = pd.read_csv(folder/"nz.csv",encoding="utf-8-sig")
    for df, start, source in [(ceh,2,"10.5285/31cc5cec-d33b-4dd6-a932-061ff947e708"),
                              (nz,1,"Manaaki Whenua 20cc42be-a0d2-439a-8b5f-ff5a7708fe3f")]:
        hits = df.iloc[:,start:].apply(pd.to_numeric).sum()
        add(hits.index[hits.gt(0)], source, "body pollen")
    dna = pd.read_csv(folder/"MothPollenFinal_prepared.csv")
    samples = pd.read_csv(folder/"moth_metadata.csv").Reference + ".blast"
    assert len(samples)==311 and samples.isin(dna.columns).all()
    present = dna.AssignmentLevel.eq("Species") & dna[samples].ge(50).any(axis=1)
    add(dna.loc[present,"Assignment"].str.replace("_"," "), "10.5281/zenodo.1408559", "body pollen DNA")
    add([r[0] for r in TABLE2 if r[3]>0], "10.1111/een.13446", "body pollen")
    spanish_tables = Document(folder/"spain.docx").tables
    for index in (0,2,3):
        for row in spanish_tables[index].rows[2:]:
            name = row.cells[0].text.strip()
            if name == "Moth species":
                break
            if int(row.cells[-1].text) > 0:
                add([" ".join(name.split()[:2])], "10.1038/s41598-023-49944-y", "body pollen")
    add(["Datura wrightii"], "10.1093/aob/mcp053", "moth-mediated pollen delivery")
    add(["Silene vulgaris"], "10.2307/2260655", "single-visit stigmatic grains")
    add(["Silene latifolia"], "10.1111/j.1365-2435.2009.01658.x", "single-visit reproductive outcome")
    add(["Platanthera chlorantha"], "Lankesteriana 17 (2017) 383-393", "pollinarium on moth")
    add(["Asclepias syriaca"], "PMC4080991", "pollinium on moth")
    return pd.DataFrame(rows,columns=["source_name","source","measurement"])


def run(inputs, signed, folder, out):
    out.mkdir(exist_ok=True,parents=True)
    s, digest, plants, _ = primary(inputs,signed)
    hosts = plants.loc[plants.all_moth].merge(s[["key","land_class"]],on="key",validate="many_to_one")
    counts = hosts.groupby("wcvp_accepted_name").land_class.nunique()
    pool = set(counts[counts.eq(2)].index)
    assert len(pool)==544
    evidence = pollen_records(folder)
    pl = pd.read_csv(folder/"gl.csv",encoding="latin1")
    pl["name"] = pl.Species_accepted_names.str.replace("_"," ").str.strip()
    sc = pd.read_excel(folder/"si.xlsx",sheet_name="Regression data")
    names = set(evidence.source_name) | set(pl.name.dropna()) | set(sc.work_species.dropna())
    names = {n for n in names if re.fullmatch(r"[A-Z][a-z]+ [a-z-]+",n)}
    crosswalk = taxonomy(inputs,names,set())
    mapped = crosswalk.loc[crosswalk.resolution.eq("resolved"),"accepted_name"].to_dict()
    mapped.update({n:n for n in pool})
    def joined(frame,col):
        frame = frame.assign(plant=frame[col].map(mapped))
        return frame[frame.plant.isin(pool)].copy()
    evidence = joined(evidence,"source_name").drop_duplicates().sort_values(["plant","source"])
    evidence.to_csv(out/"pollen_evidence.csv",index=False)
    documented = set(evidence.plant)
    pl = joined(pl,"name")
    pl = pl[np.isfinite(pl.PL_Effect_Size)].copy()
    pl[["plant","name","Author","Year","DOI","year_of_study","PL_Effect_Size","PL_Effect_Size_Type1","PL_Effect_Size_Type2"]].to_csv(out/"pollen_limitation_cases.csv",index=False)
    sc = joined(sc,"work_species")
    ledger = pd.DataFrame({"plant":sorted(pool)})
    ledger["moth_transport_documented"] = ledger.plant.isin(documented)
    ledger["sc_source_class"] = ledger.plant.map(sc.groupby("plant")["SI analyses"].agg(lambda v:"|".join(sorted(set(v))))).fillna("unknown")
    ledger["glopl_cases"] = ledger.plant.map(pl.groupby("plant").size()).fillna(0).astype(int)
    ledger["glopl_median_log_response_ratio"] = ledger.plant.map(pl.groupby("plant").PL_Effect_Size.median())
    ledger.to_csv(out/"shared_544_evidence.csv",index=False)
    components = {"transport":documented,"visit_evidence_only":pool-documented}
    for scope, group in components.items():
        c = hosts.loc[hosts.wcvp_accepted_name.isin(group)].groupby("key")[ID].nunique()
        s[scope+"_moth"] = s.key.map(c).fillna(0)
    s["shared_moth"] = s.transport_moth+s.visit_evidence_only_moth
    s["other_moth"] = s.all_moth-s.shared_moth
    assert s.other_moth.ge(0).all()
    fits = {scope:fit(s,scope,"moth",balance=True) for scope in [*components,"shared","other","all"]}
    ref=json.loads((out.parent/"shared_visitors_regional/global_composition_partition.json").read_text())
    for field in ["beta_island_minus_mainland_slope_pp_per_log1p","beta_mainland_pp_per_log1p","beta_island_pp_per_log1p"]:
        np.testing.assert_allclose(fits["transport"][field]+fits["visit_evidence_only"][field],ref["shared"][field],atol=1e-10,rtol=0)
        np.testing.assert_allclose(fits["shared"][field]+fits["other"][field],ref["all"][field],atol=1e-10,rtol=0)
    hashes = {n:hashlib.sha256((folder/n).read_bytes()).hexdigest() for n in ["rsbl20190877_si_001.xlsx","ceh.zip","nz.csv","MothPollenFinal_prepared.csv","moth_metadata.csv","spain.docx","gl.csv","si.xlsx"]}
    result = dict(target_sha256=digest,input_hashes=hashes,transport=len(documented),glopl_species=pl.plant.nunique(),
                  glopl_cases=len(pl),sc_categories=ledger.sc_source_class.value_counts().to_dict(),fits=fits)
    (out/"summary.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("inputs","signed","folder","out"):
        parser.add_argument(name,type=Path)
    run(**vars(parser.parse_args()))
