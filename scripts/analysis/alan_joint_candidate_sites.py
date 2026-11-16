from pathlib import Path
import sys
import pandas as pd

root=Path(__file__).resolve().parents[2]/"results"
table=pd.read_csv(root/"joint_night_autonomy/species.csv")
c=table[table.night_stronger & table.autonomy.eq("absent")].set_index("species")
inp=Path(sys.argv[1])
keys=pd.read_parquet(Path(sys.argv[2])/"strict_site_light.parquet")[["key"]]
sites=pd.read_pickle(inp/"sites.pkl").merge(keys,on="key",validate="1:1")
assert len(sites)==10698 and sites.region5.nunique()==78
obs=pd.read_pickle(inp/"members.pkl").merge(sites[["key","region5","year","land_class"]],on="key",validate="m:1")
obs=obs[obs.wcvp_accepted_name.isin(c.index)]
for side in ["island","mainland"]:
    c[side+"_records"]=obs.loc[obs.land_class.eq(side)].groupby("wcvp_accepted_name").key.nunique().reindex(c.index,fill_value=0)
common=obs.groupby(["wcvp_accepted_name","year","region5"]).land_class.nunique().eq(2)
c["both_region_years"]=common.groupby(level=0).sum().reindex(c.index,fill_value=0)
c.to_csv(root/"joint_night_autonomy/candidates.csv")
print(c[["night_share","island_records","mainland_records","both_region_years"]].to_string())
