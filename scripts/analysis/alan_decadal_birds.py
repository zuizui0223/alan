from alan_decadal_events import *
def birds(plant):
    rows=[]
    for offset in range(0,100000,500):
        d=get(C/f"birds_{plant}_{offset}.json","https://api.globalbioticinteractions.org/interaction",
            dict(sourceTaxon="Aves",targetTaxon=plant,format="json",limit=500,offset=offset))
        a=pd.DataFrame(d["data"],columns=d["columns"])
        rows.extend(a[a.target_taxon_name.eq(plant) & a.source_taxon_name.isin(CAT.loc[CAT.class_name.eq("Aves"),"taxon"])][["source_taxon_name","target_taxon_name","interaction_type"]].to_dict("records"))
        if len(a)<500:break
    return rows
if __name__=="__main__":
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(4) as pool: rows=sum(pool.map(birds,targets.plant.unique()),[])
    pd.DataFrame(rows).drop_duplicates().to_csv(OUT/"bird_links.csv",index=False)
