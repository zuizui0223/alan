"""Cross published reproductive function and direct autonomy on the frozen bilateral cohort."""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0,str(Path(sys.argv[1]).resolve()))
from scripts.analysis.alan_real_daynight_reproductive_PCNL_562_20261113 import experimental_daynight
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import direct_plant_traits

root=Path(__file__).resolve().parents[2]/"results"
study=Path(sys.argv[2])
d,_=experimental_daynight(study/"DvN_LiteratureDataClean_2024-09-13.csv",study/"SppTraitData_clean.csv")
d=d.set_index("species")
assert len(d)==134 and d.night_stronger.sum()==75
d["autonomy"]=direct_plant_traits(Path(sys.argv[3])).autonomous_selfing_capacity
s=pd.read_csv(root/"species_interpretation/island_sign_reproductive_evidence.csv").set_index("plant")
d=d.join(s)
d["in_shared544"]=d.island_pp.notna()
links=pd.read_csv(root/"species_interpretation/land_negative_predator_chains.csv")
d["archived_predators"]=links.groupby("plant").predator.agg(lambda x:" | ".join(sorted(set(x))))
ph=pd.read_csv(root/"long_phenology_11/coverage.csv")
d["in_regional_flowering_panel"]=d.index.isin(ph.plant)
out=root/"joint_night_autonomy";out.mkdir(exist_ok=True)
d.to_csv(out/"species.csv")
for section,g in [('all',d),('shared',d[d.in_shared544])]:
 print(section,pd.crosstab(g.autonomy.fillna('unknown'),g.night_stronger).to_string())
print(d[d.night_stronger & d.autonomy.eq('absent')][['night_share','n_studies','island_pp','mainland_pp','transport_review_status']].to_string())
print(d[d.in_shared544][['night_share','autonomy','island_pp','mainland_pp']].to_string())
