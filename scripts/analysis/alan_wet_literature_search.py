"""Catalogue searches for the focused evidence review; save identifiers, not abstracts."""
from pathlib import Path
import requests
import pandas as pd

out=Path(__file__).resolve().parents[2]/"results/all_species_forest/wet_environment"
light='(TITLE_ABS:"artificial light at night" OR TITLE_ABS:"artificial night lighting" OR TITLE_ABS:"street lighting" OR TITLE_ABS:"polarized light pollution")'
wet='(TITLE_ABS:wetland* OR TITLE_ABS:riparian OR TITLE_ABS:freshwater OR TITLE_ABS:aquatic OR TITLE_ABS:shoreline OR TITLE_ABS:pond OR TITLE_ABS:stream OR TITLE_ABS:waterfront)'
themes={"wetland_light":light+" AND "+wet,
        "plant_reproduction":light+" AND (TITLE_ABS:pollinat* OR TITLE_ABS:flowering OR TITLE_ABS:phenolog* OR TITLE_ABS:seedling*)",
        "wetland_hydrology":'(TITLE_ABS:Ellenberg AND (TITLE_ABS:"water table" OR TITLE_ABS:"water regime" OR TITLE_ABS:"soil moisture"))'}
rows=[]
queries=[]
for theme,terms in themes.items():
    query=f"({terms}) AND FIRST_PDATE:[1900-01-01 TO 2026-11-17]"
    cursor="*"
    total=0
    while True:
        res=requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={"query":query,"format":"json","pageSize":500,"cursorMark":cursor},timeout=90)
        res.raise_for_status()
        result=res.json()
        batch=result["resultList"]["result"]
        for paper in batch:
            rows.append(dict(theme=theme,source=paper["source"],id=paper["id"],year=paper.get("pubYear"),
                             doi=paper.get("doi"),title=paper.get("title")))
        total+=len(batch)
        if total>=result["hitCount"]:break
        cursor=result["nextCursorMark"]
    queries.append(dict(theme=theme,query=query,retrieved=total))
pd.DataFrame(rows).to_csv(out/"literature_search.csv",index=False)
pd.DataFrame(queries).to_csv(out/"queries.csv",index=False)
print(len(rows),pd.DataFrame(rows).drop_duplicates(["source","id"]).shape[0],[(q["theme"],q["retrieved"]) for q in queries])
