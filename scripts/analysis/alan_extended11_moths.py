from alan_extended11 import *
from collections import defaultdict
legacy=pd.read_csv(Path("../v011/old/GloBI_stage2_signed_extension.csv.gz"))
links=legacy[legacy.guild.eq("moth_family_visit")].groupby("plant_name").insect_name.agg(set)
def moth(st):
    row=coords.loc[st]
    # Restrict records to linked moth species by GBIF backbone keys, in chunks.
    spp=set.union(*(links[p] for p in targets.loc[targets.station.eq(st),"plant"]))
    keys=[]
    for name in sorted(spp):
        match=get(C/("key_"+name.replace("/","_")+".json"),"https://api.gbif.org/v1/species/match",{"name":name})
        if match.get("matchType")=="EXACT" and match.get("rank")=="SPECIES":keys.append(match.get("acceptedUsageKey",match["usageKey"]))
    for k in range(0,len(keys),40):
        for offset in range(0,100000,300):
            args=dict(taxonKey=keys[k:k+40],decimalLongitude=f"{row.longitude-.025},{row.longitude+.025}",
                decimalLatitude=f"{row.latitude-.025},{row.latitude+.025}",year="2012,2025",lifeStage="Adult",limit=300,offset=offset)
            obj=get(C/f"gbif_{st}_{k}_{offset}.json","https://api.gbif.org/v1/occurrence/search",args)
            if obj["endOfRecords"]:break
if __name__=="__main__":batch(moth,list(coords.index),5)
