"""Annual paired recording rates with a fixed local species roster and continuing survey denominator."""
from alan_wet_shore import repo,cache,frozen
import duckdb
import pandas as pd
out=repo/"results/all_species_forest/wet_annual"
db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
wet=pd.read_csv(repo/"results/all_species_forest/wet_partner_composition/species.csv")
hab=pd.read_csv(out/"habitats.csv")
water=pd.read_csv(out/"water.csv").query("status=='ok'")
db.execute("""create temp table s as select y.* from year_surveys y join water using(unit_id,y)
    semi join hab using(unit_id) qualify count(*) over(partition by site)>=5 and max(y) over(partition by site)-min(y) over(partition by site)>=20""")
db.execute("""create temp table ro as select distinct s.site,s.unit_id,h.sp,w.candidate from s join year_hits h using(site,y) join wet w using(sp)
 qualify bool_or(candidate) over(partition by site) and bool_or(not candidate) over(partition by site)""")
ro=db.table("ro").df();ro.to_csv(out/"roster.csv.gz",index=False)
raw=db.sql("""with size as (select site,count(*) filter(where candidate) nc,count(*) filter(where not candidate) nr from ro group by site),
 hits as (select h.site,h.y,sum(hits) filter(where candidate) hc,sum(hits) filter(where not candidate) hr
 from year_hits h join ro using(site,sp) group by all)
 select s.site,s.unit_id,s.y,s.block,s.visits,s.richness,s.dn,s.month_fraction_6+s.month_fraction_7+s.month_fraction_8 summer,
 100.*coalesce(hc,0)/(nc*visits) candidate_rate,100.*coalesce(hr,0)/(nr*visits) reference_rate,
 nc,nr,coalesce(hc,0)>0 and coalesce(hr,0)>0 as both_seen from s join size using(site) left join hits using(site,y)""").df()
# Require at least one co-recorded year but keep nondetections in later surveyed years.
raw["ever_both"]=raw.groupby("site").both_seen.transform("any")
raw.to_csv(out/"survey_series.csv.gz",index=False)
print(raw.groupby("unit_id").agg(yrs=("y","nunique"),series=("site","nunique"),both=("both_seen","sum"),bright=("dn","nunique")).to_string())
