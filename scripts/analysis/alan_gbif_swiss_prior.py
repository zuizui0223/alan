"""Compare Swiss grid results with MOIS and available reproductive measurements."""
from pathlib import Path
import sys
import pandas as pd

out=Path(sys.argv[1])
res=pd.read_csv(out/"gbif_species.csv")
prior=pd.read_csv(out.parent/"register_long_light/species.csv")
cols=["species","time_pp_per_decade","q_time","light_pp_per_10dn","q_light","autonomous_selfing_capacity","self_incompatibility","night_share","glopl_median_log_response_ratio"]
prior=prior[cols].dropna(subset=["species"])
assert prior.species.is_unique
d=res.merge(prior.rename(columns={"q_time":"mois_q_time","q_light":"mois_q_light"}),how="left",on="species",validate="1:1")
d.to_csv(out/"gbif_mois_evidence.csv",index=False)
prev=prior[prior.q_time.lt(.05)&prior.time_pp_per_decade.lt(0)]
comparison=prev[["species","time_pp_per_decade"]].merge(res,how="left",on="species")
comparison.to_csv(out/"previous_three.csv",index=False)
print(comparison.to_string(index=False))
print("tested",d.q_time.count(),d.q_light.count(),"declines",sum(d.q_time.lt(.05)&d.delta_pp.lt(0)),"negative light",sum(d.q_light.lt(.05)&d.pp_per_10dn.lt(0)))
for title,g in [("tested",d[d.q_time.notna()]),("decliners",d[d.q_time.lt(.05)&d.delta_pp.lt(0)])]:
 print(title,len(g),g.night_share.notna().sum(),g.glopl_median_log_response_ratio.notna().sum(),g.autonomous_selfing_capacity.value_counts().to_dict())
