from alan_extended11 import coords, targets, phases, names, get, batch,F
from alan_decadal_light import C,OUT
import pandas as pd
import numpy as np
import requests
CAT=pd.read_json(F/"npn/species.json")
CAT["taxon"]=CAT.genus+" "+CAT.species
LEG=pd.read_csv("../v011/old/GloBI_stage2_signed_extension.csv.gz")
LEG=LEG[LEG.plant_name.isin(targets.plant)]
IDS=set(names[names.accepted_name.isin(targets.plant)].species_id)
LEP=set(CAT.loc[CAT.taxon.isin(LEG.insect_name),"species_id"])
BIRDS=set(CAT.loc[CAT.class_name.eq("Aves"),"species_id"])
PH=pd.read_json(F/"npn/getPhenophases.json")
def fetch(job):
    year,kind=job
    stations=list(coords.index)
    params=dict(start_date=f"{year}-01-01",end_date=f"{year}-12-31",request_src="Historical phenology research")
    ids=IDS if kind=="plant" else LEP if kind=="insect" else BIRDS
    phen=phases if kind=="plant" else PH.loc[PH.phenophase_name.str.strip().str.lower().isin(["adults","winged adults","adults feeding","individuals","feeding"]),"phenophase_id"].tolist()
    for key,group in [("station_id",stations),("species_id",sorted(ids)),("phenophase_id",phen)]:
        params.update({f"{key}[{i}]":int(x) for i,x in enumerate(group,1)})
    d=get(C/f"npn_{year}_{kind}.json","https://services.usanpn.org/npn_portal/observations/getObservations.json",params,True)
    assert isinstance(d,list) and all(x["species_id"] in ids and x["site_id"] in stations for x in d)
def biomass(job):
    st,y=job
    r=coords.loc[st]
    item=f"chloris_biomass_50km_{y}"
    d=get(C/f"woody_{st}_{y}.json",f"https://planetarycomputer.microsoft.com/api/data/v1/item/point/{r.longitude},{r.latitude}",
          dict(collection="chloris-biomass",item=item,assets="biomass"))
    assert d["band_names"]==["biomass_b1"]
if __name__=="__main__":
    batch(fetch,[(y,k) for y in range(2008,2025) for k in ["plant","insect","bird"]],6)
    batch(biomass,[(st,y) for st in coords.index for y in range(2003,2020)],12)
