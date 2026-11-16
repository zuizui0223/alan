"""Review European riparian floral contacts for the frozen wet-soil cohort.

Download Data/3_Final_data/Interaction_data.csv.gz from JoseBSL/EuPPollNet
tag v1.3.0 to work/gbif_annual/euppollnet_v130.csv.gz before running.
"""
from pathlib import Path
from zipfile import ZipFile

import duckdb
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from alan_rf_decline_traits import annotate, compare

repo = Path(__file__).resolve().parents[2]
work = repo.parent
results = repo / "results"
root = results / "all_species_forest"
out = root / "wet_pollination"
out.mkdir(exist_ok=True)
wet = pd.read_csv(root / "decline_traits/moisture_evidence.csv.gz").query("indicator > 7").set_index("sp")
assert (len(wet), wet.candidate.sum()) == (363, 58)
traits, _, wet["accepted"] = annotate(work, results, wet)
cols = ["Study_id", "Network_id", "Date", "Latitude", "Longitude", "Interaction",
        "EuPPollNet_habitat", "Plant_accepted_name", "Plant_unsure_id", "Plant_rank",
        "Pollinator_accepted_name", "Pollinator_unsure_id"]
net = pd.read_csv(work / "gbif_annual/euppollnet_v130.csv.gz", encoding="latin1", usecols=cols)
net = net.loc[net.Plant_unsure_id.eq("No") & net.Plant_rank.eq("SPECIES") & net.Interaction.gt(0)].copy()
with ZipFile(work / "v011/wcvp.zip") as archive:
    names = pd.read_csv(archive.open("wcvp_names.csv"), sep="|", dtype=str,
                        usecols=["taxon_name","taxon_rank","taxon_status","plant_name_id","accepted_plant_name_id"])
names = names[names.taxon_rank.eq("Species")]
accepted = names[names.taxon_status.eq("Accepted")].set_index("plant_name_id").taxon_name
aliases = names[names.taxon_name.isin(net.Plant_accepted_name)].copy()
aliases["target"] = aliases.accepted_plant_name_id.where(~aliases.taxon_status.eq("Accepted"), aliases.plant_name_id)
aliases["valid"] = aliases.taxon_status.isin(["Accepted","Synonym","Heterotypic Synonym","Homotypic Synonym"]) & aliases.target.isin(accepted.index)
groups = aliases.groupby("taxon_name")
mapping = groups.target.first().where(groups.valid.all() & groups.target.nunique().eq(1)).map(accepted)
net["accepted"] = net.Plant_accepted_name.map(mapping)
net = net[net.accepted.isin(wet.accepted.dropna())].copy()
net["riparian"] = net.EuPPollNet_habitat.isin(["Riparian vegetation", "Riparian"])
net["dronefly"] = net.Pollinator_accepted_name.eq("Eristalis tenax") & net.Pollinator_unsure_id.eq("No")
for label, subset in [("any",net),("riparian",net[net.riparian]),("dronefly",net[net.dronefly]),
                      ("riparian_dronefly",net[net.riparian & net.dronefly])]:
    wet[label+"_studies"] = wet.accepted.map(subset.groupby("accepted").Study_id.nunique()).fillna(0).astype(int)
selected = ["moth_visit_documented","moth_transport_reviewed","night_exclusion_loss","pollen_limitation"]
wet = wet.join(traits[selected])
wet.to_csv(out / "species.csv")
rows = []
for col in selected + ["riparian_studies","dronefly_studies","riparian_dronefly_studies"]:
    v = wet[col]
    if col.endswith("_studies"):
        coverage = wet.riparian_studies if col == "riparian_dronefly_studies" else wet.any_studies
        v = v.gt(0).astype(float).where(coverage.gt(0))
    row = compare(wet,v,["site_region","island_band","site_band"],np.random.default_rng(8102))
    rows.append(dict(trait=col, controls_known=v[~wet.candidate].count(), **row))
pd.DataFrame(rows).to_csv(out / "contrasts.csv",index=False)
# Distance to a water-analysis plant survey in the same calendar year, regardless of species.
db = duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
xy = db.sql("select distinct unit_id, latitude, longitude from geo").df()
survey = pd.read_csv(root/"wet_water/satellite_survey_years.csv.gz").merge(xy,validate="m:1")
survey = survey[survey.unit_id.isin(pd.read_csv(root/"wet_water/coordinate_trends.csv.gz").unit_id)]
sites = net[net.riparian].copy()
sites["year"] = pd.to_numeric(sites.Date.str[:4],errors="coerce")
for year,g in sites.groupby("year"):
    ref = survey[survey.y.eq(year)]
    if len(ref):
        tree = BallTree(np.radians(ref[["latitude","longitude"]]),metric="haversine")
        dist,_ = tree.query(np.radians(g[["Latitude","Longitude"]]),k=1)
        sites.loc[g.index,"nearest_survey_km"] = dist[:,0]*6371.009
sites.to_csv(out/"riparian_contacts.csv.gz",index=False)
