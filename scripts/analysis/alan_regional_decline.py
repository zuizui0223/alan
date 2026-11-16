"""Prioritize exploratory regional signals using time trends and nearby non-brightening sites."""
from pathlib import Path
import sys
import duckdb
import pandas as pd
import numpy as np

def run(cache,out):
    d=pd.read_csv(out/"all_regions.csv.gz")
    db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
    g=db.sql("""with xy as (select distinct unit_id,longitude,latitude,land_class,land_id from geo)
       select t.*, case when t.land_class='island' then 'I_'||xy.land_id else
       'M_'||cast(floor((longitude+180)/5) as int)||'_'||cast(floor((latitude+90)/5) as int) end as region
       from trajectories t join xy using(unit_id,land_class)""").df()
    assert not g.duplicated(["unit_id","sp"]).any()
    # A near-zero fitted slope from floating point alone is assigned to the reference.
    g["bright"]=g.delta>1e-8
    blocks=g.groupby(["sp","region","block","bright"],as_index=False).agg(rate=("change","mean"),dn=("delta","mean"),coords=("unit_id","size"))
    paired=blocks.pivot(index=["sp","region","block"],columns="bright",values=["rate","dn","coords"]).dropna()
    x=pd.DataFrame(dict(bright_rate=paired["rate"][True],ref_rate=paired["rate"][False],
                        bright_dn=paired["dn"][True],ref_dn=paired["dn"][False],
                        bright_n=paired["coords"][True],ref_n=paired["coords"][False])).reset_index()
    x["difference"]=x.bright_rate-x.ref_rate
    x["joint"]=(x.bright_rate<0)&(x.difference<0)
    sm=x.groupby(["sp","region"],as_index=False).agg(paired_blocks=("block","size"),joint_blocks=("joint","sum"),
         bright_n=("bright_n","sum"),ref_n=("ref_n","sum"),bright_rate=("bright_rate","mean"),ref_rate=("ref_rate","mean"),
         difference=("difference","mean"),bright_dn=("bright_dn","mean"),ref_dn=("ref_dn","mean"))
    y=d[d.label.eq("neg")].merge(sm,how="left",validate="1:1")
    y["priority"]=(y.paired_blocks>=3)&(y.joint_blocks>y.paired_blocks/2)&(y.bright_rate<0)&(y.difference<0)
    y.to_csv(out/"negative_temporal_check.csv",index=False)
    traits = pd.read_csv(out/"species_traits.csv")
    shortlist = y.loc[y.priority].merge(traits, on=["sp", "species"], validate="m:1")
    assert len(shortlist) == y.priority.sum()
    shortlist.sort_values(["joint_blocks", "paired_blocks", "species"], ascending=[False, False, True]).to_csv(
        out/"priority_species_regions.csv", index=False)
    print(y.priority.sum(),y.loc[y.priority,"sp"].nunique(),y[y.priority].region.nunique())
    print(y[y.priority].sort_values(["joint_blocks","paired_blocks"],ascending=False)[["species","region","bright_rate","ref_rate","joint_blocks","paired_blocks"]].head(15).round(2).to_string(index=False))
if __name__=="__main__":run(*map(Path,sys.argv[1:3]))
