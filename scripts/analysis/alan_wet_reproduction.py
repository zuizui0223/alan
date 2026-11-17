"""Reconcile wet-cohort mating, vector and moth-spider evidence, retaining missing states."""
from pathlib import Path
import io
import subprocess

import duckdb
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from alan_rf_decline_traits import compare

repo=Path(__file__).resolve().parents[2]
cache=repo.parent/"gbif_annual"
out=repo/"results/all_species_forest/wet_reproduction"
out.mkdir(exist_ok=True)
def frozen(file):
    b=subprocess.check_output(["git","show","b780370:results/all_species_forest/wet_network/"+file],cwd=repo)
    return pd.read_csv(io.BytesIO(b),compression="gzip" if file.endswith(".gz") else None)
w=frozen("species.csv").set_index("sp")
a=pd.read_csv(cache/"pladias_wet_accounts.csv").set_index("sp")
assert set(a.index)==set(w.index) and w.candidate.sum()==58 and len(w)==363
w=w.join(a)
for col in ["Pollination syndrome","Generative reproduction type","Reproduction type","Flower colour","Flower symmetry"]:
    w[col]=w[col].str.replace(r"\s+"," ",regex=True).str.strip().replace(r"^\s*[-–—?]?\s*$",np.nan,regex=True)
def has(text,pattern):
    return text.str.contains(pattern,case=False,regex=True).where(text.notna()).astype(float)
v=w["Pollination syndrome"]
features={}
for label,pat in [("insect_vector","insect-pollination"),("wind_vector","wind-pollination"),("water_vector","water-pollination"),
                  ("self_transfer","selfing|cleistogamy|geitonogamy")]:
    features[label]=has(v,pat)
breed=w["Generative reproduction type"]
features["mixed_mating"]=has(breed,"mixed mating")
features["autogamy_or_facultative"]=has(breed,r"\bautogamy\b")
features["allogamy_or_facultative"]=has(breed,r"\ballogamy\b")
features["vegetative_capacity"]=has(w["Reproduction type"],"vegetativ")
# Local direct-trait ledger: summarize restricted traits only in aggregate.
ledger=pd.read_csv(repo.parent/"direct_trait_replay/direct_species_trait_ledger.csv.gz",low_memory=False)
ledger=ledger[ledger.resolution_status.eq("resolved") & ledger.quality.isin(["high","medium"]) & ledger.multistate.eq(False)]
traits=ledger.pivot(index="accepted_species",columns="trait_name",values="normalized_value").reindex(w.accepted).set_axis(w.index)
for label,col,yes,no in [("self_incompatible","self_incompatibility","SI","SC"),
                         ("autonomous_selfing","autonomous_selfing_capacity","autonomous","absent")]:
    s=traits[col]
    features[label]=s.eq(yes).where(s.isin([yes,no])).astype(float)
features["white_tubular_proxy"]=(traits.flower_primary_color.eq("white") & traits.floral_form.isin(["tubular","salverform","funnel_trumpet"])).where(traits.flower_primary_color.notna() & traits.floral_form.notna()).astype(float)
features["moth_visit_non_inat"]=w.moth_visit_documented
features["moth_transport_reviewed"]=w.moth_transport_reviewed

# Require named non-butterfly Lepidoptera with an unambiguous, identified family.
fl=frozen("globi_flower_sources.csv.gz")
diet=frozen("diet_sources.csv.gz")
db=duckdb.connect()
db.read_parquet(str(cache/"globi_v011.parquet")).create_view("src")
def tax(side):
    return f"""select distinct coalesce(nullif({side}TaxonSpeciesName,''),{side}TaxonName) as visitor,
      {side}TaxonFamilyName as family from src where {side}TaxonOrderName='Lepidoptera'
      and regexp_full_match({side}TaxonFamilyName,'[A-Z][a-z]+idae')"""
taxa=db.sql(tax("source")+" union "+tax("target")).df()
fam=taxa.groupby("visitor").family
families=fam.first().where(fam.nunique().eq(1)).dropna()
butterfly={"Papilionidae","Nymphalidae","Pieridae","Lycaenidae","Riodinidae","Hesperiidae","Hedylidae"}
fl=fl[fl.visitor.isin(families[~families.isin(butterfly)].index)]
spider=diet[diet.guild.eq("spiders") & diet.prey_stage.str.lower().isin(["adult","adults","adult(s)","imago"])
           & diet.dataset.notna() & ~diet.dataset.str.lower().str.contains("trophich",na=False)]
chain=fl.merge(spider,on="visitor",suffixes=("_flower","_diet"))
chain.to_csv(out/"moth_spider_sources.csv.gz",index=False)
wetflowers=frozen("flower_sources.csv.gz")
wetflowers=wetflowers[wetflowers.wet_site & wetflowers.visitor.isin(families[~families.isin(butterfly)].index)]
wet_chain=wetflowers.merge(spider,on="visitor")
wet_chain.to_csv(out/"habitat_labelled_moth_spiders.csv",index=False)
for label,group in [("any",chain),("non_inat_flower",chain[chain.dataset_flower.notna() & ~chain.dataset_flower.str.lower().str.contains("inaturalist",na=False)]),
                    ("wet_flower",wet_chain)]:
    w["moth_spider_"+label]=w.accepted.isin(group.plant)
features["moth_spider_among_named_visitors"]=w.moth_spider_any.astype(float).where(w.named_flower)

rows=[]
for name,value in features.items():
    result=compare(w,value,["site_region","island_band","site_band"],np.random.default_rng(7125))
    rows.append(dict(trait=name,candidate_yes=value[w.candidate].sum(),candidate_known=value[w.candidate].count(),
                     reference_yes=value[~w.candidate].sum(),reference_known=value[~w.candidate].count(),**result))
summary=pd.DataFrame(rows)
usable=summary.p.notna()
summary.loc[usable,"holm_p"]=multipletests(summary.loc[usable,"p"],method="holm")[1]
summary.to_csv(out/"contrasts.csv",index=False)
spectrum=pd.read_excel(cache/"pladias_320.xlsx").groupby("Taxon")["Value [enum]"].agg(lambda x:", ".join(sorted(set(x.dropna())-{"unknown"})))
w["documented_visitor_groups"]=w.account_name.map(spectrum)
# Walton et al. S3 gives counts by moth family for each named pollen source.
walton=pd.read_excel(cache/"walton_rsbl20190877_si_001.xlsx",sheet_name="Plant-Pollinator Table (S3)",header=3)
name=walton.iloc[:,0].str.strip().str.rstrip("*").replace({"Iris psuedacorus":"Iris pseudacorus"})
families_s3=["Noctuidae","Geometridae","Erebidae","Crambidae","Hepliaidae","Lasiocampidae","Notodontidae"]
counts=walton[families_s3].sum(axis=1,min_count=1)
evidence=pd.Series(counts.to_numpy(),index=name).dropna()
assert evidence.index.is_unique
w["walton_s3_moths"]=w.accepted.map(evidence).combine_first(w.species.map(evidence))
w["documented_transport"]=w.moth_transport_reviewed.eq(1) | w.walton_s3_moths.gt(0)
w[["species","accepted","candidate","indicator","url","account_name","status","Pollination syndrome","Generative reproduction type","Reproduction type","Flower colour","Flower symmetry","documented_visitor_groups","moth_visit_documented","moth_transport_reviewed","night_exclusion_loss","pollen_limitation","walton_s3_moths","documented_transport","moth_spider_any","moth_spider_non_inat_flower","moth_spider_wet_flower"]].to_csv(out/"species.csv")
