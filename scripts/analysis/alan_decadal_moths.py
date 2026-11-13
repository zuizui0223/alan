from alan_extended11_moths import moth, coords, batch
from alan_decadal_light import C,OUT
from pathlib import Path
import json
import pandas as pd
from alan_decadal_events import LEG, targets
if __name__=="__main__":
    batch(lambda st:moth(st,1992,2011,C),list(coords.index),5)
    old=pd.read_csv(OUT.parent/"long_phenology_11/moth_dates.csv")
    rows=[]
    for t in targets.itertuples():
        names=set(LEG.loc[LEG.guild.eq("moth_family_visit")&LEG.plant_name.eq(t.plant),"insect_name"])
        for path in C.glob(f"gbif_{t.station}_*.json"):
            for r in json.loads(path.read_bytes())['results']:
                if r.get('species') in names and r.get('day') and r.get('lifeStage')=='Adult':
                    rows.append(dict(station=t.station,plant=t.plant,land=t.land_class,moth=r['species'],date=r['eventDate'][:10],gbif_key=r['key']))
    d=pd.concat([old[old.date.str[:4].astype(int)<=2024],pd.DataFrame(rows)]).drop_duplicates()
    d.to_csv(OUT/'gbif_moth_dates.csv',index=False)
