"""Join the publisher's GB trend ensembles, retaining unknown traits and unreliable trends."""
import sys
import csv
from pathlib import Path
from zipfile import ZipFile
import pandas as pd

cache, out, pr85, ledger = map(Path, sys.argv[1:5])
sys.path.insert(0,str(pr85.resolve()))
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import direct_plant_traits

trends = cache / "PA2020_HectadTrendsData"
res = pd.read_csv(out/"atlas_light.csv").set_index("id")
mapped = pd.read_csv(trends/"aggregateMappings.csv")
def read_one(archive,path):
    return next(csv.DictReader(archive.read(path).decode().splitlines()))
with ZipFile(trends/"hectad-counts.zip") as counts:
    for period in ["long","short"]:
        exclude = set(mapped.loc[mapped.analysisType.eq(period),"mapped.ddb.id"].str.replace(".","_",regex=False))
        with ZipFile(trends/f"{period}.zip") as z:
            for file in z.namelist():
                if not file.startswith(f"{period}/trends-summaries/Britain/") or not file.endswith(".csv"):
                    continue
                key = Path(file).stem
                if key not in res.index or key in exclude:
                    continue
                n = int(read_one(counts,f"hectad-counts/{key}.csv")["Gb"+period.title()])
                if n<16:
                    continue
                record = read_one(z,file)
                res.loc[key,period+"_decline_of_100"] = int(record["declineMod"])+int(record["declineStrong"])
res = res.join(direct_plant_traits(ledger),on="species")
experiments = pd.read_csv(out.parent/"joint_night_autonomy/species.csv").set_index("species")
res = res.join(experiments[["night_share","night_stronger"]],on="species")
shared = pd.read_csv(out.parent/"species_interpretation/island_sign_reproductive_evidence.csv").set_index("plant")
res = res.join(shared[["transport_review_status","glopl_median_log_response_ratio","island_pp","mainland_pp"]],on="species")
res.to_csv(out/"species_evidence.csv")
decline = res.short_decline_of_100.ge(95)
both = decline & res.long_decline_of_100.ge(95)
light = res.q.lt(.05)&res.pp_per_10dn.lt(0)
print(res[light][["species","short_decline_of_100","long_decline_of_100","autonomous_selfing_capacity","self_incompatibility","night_share"]].to_string())
for flag in [res.short_decline_of_100.notna(),decline,both,decline & res.q.notna(),decline & light]:
 d=res[flag]
 print(len(d), d.autonomous_selfing_capacity.value_counts().to_dict(),d.night_stronger.value_counts().to_dict())
print(res[both & res.p.notna()].nsmallest(12,"pp_per_observed_iqr")[["species","q","pp_per_observed_iqr"]].to_string())
eligible = res[res.short_decline_of_100.notna()].copy()
eligible["declining"] = eligible.short_decline_of_100.ge(95)
comparison = eligible.groupby("declining").agg(
    species=("species","size"),
    light_fits=("q","count"),
    negative_light_estimate=("pp_per_10dn", lambda s:s.lt(0).sum()),
    median_pp_per_iqr=("pp_per_observed_iqr","median"),
    measured_absent=("autonomous_selfing_capacity", lambda s:s.eq("absent").sum()),
    measured_autonomous=("autonomous_selfing_capacity", lambda s:s.eq("autonomous").sum()),
    diel_experiment=("night_share","count"),
    night_higher=("night_share",lambda s:s.gt(.5).sum()),
    moth_transport_documented_in_shared_cohort=("transport_review_status",lambda s:s.eq("documented").sum()),
    pollen_limitation_available_in_shared_cohort=("glopl_median_log_response_ratio","count"),
)
comparison.to_csv(out/"trait_comparison.csv")
print(comparison.to_string())
