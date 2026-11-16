"""Rebuild species-matched wetland interaction inventories from frozen primary data.

Inputs in work/gbif_annual: globi_v011.parquet (Zenodo 22691479), dopi.zip
(15224385), eupolltrait.zip (18032357), euppollnet_v130.csv.gz (15183272).
Run with the previous wet_pollination species table and WCVP v011 archive.
"""
from pathlib import Path
from zipfile import ZipFile
from io import BytesIO
import hashlib
import sys
from tempfile import TemporaryDirectory

import duckdb
import numpy as np
import pandas as pd

from alan_rf_decline_traits import compare

repo = Path(__file__).resolve().parents[2]
work = repo.parent
cache = work / "gbif_annual"
out = repo / "results/all_species_forest/wet_network"
site_year_mode = "--site-years" in sys.argv
if site_year_mode:
    temporary = TemporaryDirectory()
    output_root = out.parent
    out = Path(temporary.name)
out.mkdir(exist_ok=True)
plants = pd.read_csv(repo / "results/all_species_forest" / "wet_pollination/species.csv").set_index("sp")
assert len(plants) == 363 and plants.candidate.sum() == 58
with ZipFile(work / "v011/wcvp.zip") as z:
    fields = ["taxon_name", "taxon_rank", "taxon_status", "plant_name_id", "accepted_plant_name_id"]
    names = pd.read_csv(z.open("wcvp_names.csv"), sep="|", dtype=str, usecols=fields)
names = names[names.taxon_rank.eq("Species")].copy()
valid = names.taxon_status.eq("Accepted")
accepted = names[valid].set_index("plant_name_id").taxon_name
names["resolved"] = names.accepted_plant_name_id.where(~valid,names.plant_name_id).map(accepted)
possible = names.loc[names.resolved.isin(plants.accepted.dropna()), "taxon_name"]
names = names[names.taxon_name.isin(possible)].copy()
names["valid"] = names.taxon_status.isin(["Accepted","Synonym","Heterotypic Synonym","Homotypic Synonym"]) & names.resolved.notna()
g = names.groupby("taxon_name")
mapping = g.resolved.first().where(g.valid.all() & g.resolved.nunique().eq(1)).dropna()
aliases = mapping.rename_axis("alias").reset_index(name="plant")

with ZipFile(cache / "eupolltrait.zip") as z:
    sheet = pd.read_excel(BytesIO(z.read(next(k for k in z.namelist() if k.endswith("traits_raw_bees_hoverflies.xlsx")))),header=None)
assert list(sheet.iloc[2,31:33]) == ["on_water_plants","submerged_sediment"]
aquatic = set(sheet.loc[sheet[1].eq("Syrphidae") & (sheet[31].eq("yes") | sheet[32].eq("yes")),4])
with ZipFile(cache / "dopi.zip") as outer:
    with ZipFile(BytesIO(outer.read("DoPI_dataset/DoPI_dataset.zip"))) as inner:
        d = pd.read_csv(inner.open("DoPI_dataset.csv"),encoding="latin1",low_memory=False)
d = d[d.Pollination.isin([1,2,3])].copy()
d["plant"] = d["Plant Species"].map(mapping)
d["visitor"] = d["Pollinator Species"]
d["habitat"] = d.Habitat.fillna("")
d["wet_site"] = d.habitat.isin(["Riparian Zone","Fen, Marsh and Swamp","Saltmarsh"])
d["strong"] = d.Pollination.eq(1)
d["study"] = d.Title
d["archive"] = "DoPI"
e = pd.read_csv(cache/"euppollnet_v130.csv.gz",encoding="latin1",low_memory=False)
e = e[e.Plant_unsure_id.eq("No") & e.Plant_rank.eq("SPECIES") & e.Interaction.gt(0)].copy()
e["plant"] = e.Plant_accepted_name.map(mapping)
e["visitor"] = e.Pollinator_accepted_name.where(e.Pollinator_unsure_id.eq("No") & e.Pollinator_rank.eq("SPECIES"))
e["habitat"] = e.EuPPollNet_habitat
e["wet_site"] = e.habitat.isin(["Riparian vegetation","Riparian"])
e["strong"] = False
e["study"] = e.Study_id
e["archive"] = "EuPPollNet"
cols = ["archive","study","habitat","plant","visitor","wet_site","strong"]
flowers = pd.concat([d[cols],e[cols]]).dropna(subset="plant").drop_duplicates()
flowers["aquatic_hoverfly"] = flowers.visitor.isin(aquatic)
flowers.to_csv(out/"flower_sources.csv.gz",index=False)
for label, subset in {"flower_db":flowers,"wet_flower":flowers[flowers.wet_site],
                      "wet_strong":flowers[flowers.wet_site & flowers.strong],
                      "aquatic_hoverfly":flowers[flowers.aquatic_hoverfly],
                      "wet_aquatic_hoverfly":flowers[flowers.wet_site & flowers.aquatic_hoverfly]}.items():
    plants[label] = plants.accepted.isin(subset.plant)

db = duckdb.connect()
db.execute("set threads=4")
path = cache / "globi_v011.parquet"
h = hashlib.md5()
with path.open("rb") as src:
    for block in iter(lambda:src.read(8*1024**2),b""):
        h.update(block)
assert h.hexdigest()=="8212776b7ba0c4fa8de26a00b5ce8595"
db.read_parquet(str(path)).create_view("globi")
def tax(side):
    return f"coalesce(nullif({side}TaxonSpeciesName,''),{side}TaxonName)"
def select_direct(side, partner):
    return f"""select plant, '{side}' as plant_side, interactionTypeName as relation,
      {tax(partner)} as partner, {partner}TaxonClassName as partner_class,
      {partner}TaxonOrderName as partner_order, sourceNamespace as dataset,
      referenceDoi as doi, decimalLatitude as lat, decimalLongitude as lon,eventDate as event_date
      from globi join aliases on {tax(side)}=alias
      where coalesce({side}TaxonKingdomName,'') in ('','Plantae','Viridiplantae')"""
db.execute("create temp table direct as "+select_direct("source","target")+" union all "+select_direct("target","source"))
direct = db.sql("""select plant,plant_side,relation,partner,partner_class,partner_order,
    count(*) as source_rows,count(distinct dataset) as datasets from direct group by all""").df()
direct.to_csv(out/"direct_inventory.csv.gz",index=False)
direct_summary = direct.assign(candidate=direct.plant.isin(plants.loc[plants.candidate,"accepted"]))
direct_summary.groupby(["candidate","plant_side","relation"]).agg(
    plants=("plant","nunique"),named_partners=("partner","nunique"),
    rows=("source_rows","sum")).reset_index().to_csv(out/"direct_relations.csv",index=False)
db.execute("""create temp view floral as select * from direct where partner_class='Insecta' and
  ((plant_side='target' and relation in ('visitsFlowersOf','pollinates')) or
   (plant_side='source' and relation in ('flowersVisitedBy','pollinatedBy')))
  and regexp_full_match(partner,'[A-Z][a-zA-Z-]+ [a-z][a-z-]+')""")
# Restrict external records to confidently named insect visitors for the trophic join.
# DoPI covers flower-visiting insects but retains other visitors as well; require an
# insect classification in GloBI on the prey side for every predation link.
ext = flowers[flowers.visitor.str.fullmatch(r"[A-Z][a-zA-Z-]+ [a-z][a-z-]+",na=False)]
db.execute("""create temp table links as select distinct plant,partner as visitor,'GloBI' as archive from floral
    union select distinct plant,visitor,archive from ext""")
db.execute("create temp table visitors as select distinct visitor from links")
def select_diet(pred, prey, relations):
    return f"""select visitor,{tax(pred)} as predator,{pred}TaxonOrderName as pred_order,
    {pred}TaxonClassName as pred_class, {prey}LifeStageName as prey_stage,
    interactionTypeName as relation,sourceNamespace as dataset,referenceDoi as doi,
    decimalLatitude as lat,decimalLongitude as lon,eventDate as event_date
    from globi join visitors on {tax(prey)}=visitor
    where {prey}TaxonClassName='Insecta' and interactionTypeName in ({relations})
    and {pred}TaxonKingdomName='Animalia' and regexp_full_match({tax(pred)},'[A-Z][a-zA-Z-]+ [a-z][a-z-]+')"""
db.execute("create temp table diet as "+select_diet("source","target","'eats','preysOn'")+" union all "+select_diet("target","source","'eatenBy','preyedUponBy'"))
db.execute("""create temp view pred as select *,case when pred_order='Chiroptera' then 'bats'
    when pred_class='Amphibia' then 'amphibians' when pred_class='Reptilia' or pred_order in ('Squamata','Testudines','Crocodylia','Rhynchocephalia') then 'reptiles'
    when pred_order='Odonata' then 'dragonflies_damselflies' when pred_order='Araneae' then 'spiders'
    when pred_order='Mantodea' then 'mantises' when pred_class='Aves' then 'birds'
    when pred_class='Insecta' then 'other_insects' else 'other_animals' end as guild from diet""")
diet = db.sql("""select visitor,predator,guild,prey_stage,relation,dataset,doi,count(*) as source_rows from pred group by all""").df()
diet.to_csv(out/"diet_sources.csv.gz",index=False)
db.sql("select * from links order by all").df().to_csv(out/"flower_links.csv.gz",index=False)
db.sql("select distinct plant,partner as visitor,dataset,doi,relation from floral").df().to_csv(out/"globi_flower_sources.csv.gz",index=False)
direct_set = set(direct.plant)
plants["globi_direct"] = plants.accepted.isin(direct_set)
plants["globi_flower"] = plants.accepted.isin(db.sql("select distinct plant from floral").df().plant)
all_flower = db.sql("select distinct plant from links").df()
plants["named_flower"] = plants.accepted.isin(all_flower.plant)
con = db.sql("""select distinct plant,guild from links join pred using(visitor)""").df()
for guild,group in con.groupby("guild"):
    plants[guild] = plants.accepted.isin(group.plant)
wetpairs = ext[ext.wet_site][["plant","visitor"]].drop_duplicates()
chains = db.sql("""select w.plant,w.visitor,p.guild,count(distinct p.predator) as predator_species,
  count(distinct case when lower(p.prey_stage) in ('adult','adults','adult(s)','imago') then p.predator end) as adult_prey_predators,
  count(distinct case when lower(p.prey_stage) in ('adult','adults','adult(s)','imago')
    and p.dataset is not null and not contains(lower(p.dataset),'trophich') then p.predator end) as adult_other_sources
  from wetpairs w join pred p using(visitor) group by all""").df()
chains.to_csv(out/"wet_flower_diet.csv.gz",index=False)
wet_detail = chains.merge(wetpairs,validate="m:1")
wet_detail["aquatic_hoverfly"] = wet_detail.visitor.isin(aquatic)
wet_detail["candidate"] = wet_detail.plant.isin(plants.loc[plants.candidate,"accepted"])
evidence = []
for (candidate,guild),g in wet_detail.groupby(["candidate","guild"]):
    multi = g.groupby("plant").visitor.nunique()
    evidence.append(dict(candidate=candidate,guild=guild,plants=g.plant.nunique(),visitors=g.visitor.nunique(),
                         plants_two_visitors=int(multi.ge(2).sum()),
                         plants_without_honeybee=g.loc[g.visitor.ne("Apis mellifera"),"plant"].nunique(),
                         plants_adult_prey=g.loc[g.adult_prey_predators.gt(0),"plant"].nunique(),
                         plants_adult_outside_trophich=g.loc[g.adult_other_sources.gt(0),"plant"].nunique(),
                         plants_aquatic_hoverfly=g.loc[g.aquatic_hoverfly,"plant"].nunique()))
pd.DataFrame(evidence).to_csv(out/"wet_guild_evidence.csv",index=False)
for guild,group in chains.groupby("guild"):
    plants["wet_"+guild] = plants.accepted.isin(group.plant)
# Preserve exact dated visits, including the two habitat-labelled archives.
d["event_date"] = pd.to_datetime(dict(year=d.Year,month=d.Month,day=d.Date),errors="coerce").dt.strftime("%Y-%m-%d")
e["event_date"] = e.Date
keys = ["plant","visitor","Latitude","Longitude","event_date","archive","wet_site"]
events = pd.concat([d[keys],e[keys]]).dropna(subset=["plant","visitor","Latitude","Longitude","event_date"]).drop_duplicates()
db.register("events",events)
db.execute("""create temp view dated_visits as select plant,partner,lat,lon,event_date,dataset,
    null::boolean as wet_site from floral union all select plant,visitor,Latitude,Longitude,event_date,
    archive,wet_site from events""")
if site_year_mode:
    from alan_wet_site_years import run
    run(db,d,plants,work,output_root/"wet_site_years")
    sys.exit(0)
# Exact dated/located records from GloBI only. Identical georeferences can be shared study centroids.
same = db.sql("""select distinct f.plant,f.partner as visitor,p.predator,p.guild,f.lat,f.lon,
    substr(f.event_date,1,10) as day,f.dataset as floral_dataset,p.dataset as diet_dataset,p.prey_stage,f.wet_site
    from dated_visits f join pred p on f.partner=p.visitor and try_cast(f.lat as double)=try_cast(p.lat as double)
    and try_cast(f.lon as double)=try_cast(p.lon as double) and substr(f.event_date,1,10)=substr(p.event_date,1,10)
    where regexp_matches(f.event_date,'^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )') and
    regexp_matches(p.event_date,'^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )') and
    abs(try_cast(f.lat as double))<=90 and abs(try_cast(f.lon as double))<=180 and
    (try_cast(f.lat as double)<>0 or try_cast(f.lon as double)<>0)""").df()
same.to_csv(out/"same_coordinate_day.csv",index=False)
wet_degree = wetpairs.groupby("plant").visitor.nunique()
plants["wet_visitor_richness"] = plants.accepted.map(wet_degree).fillna(0).astype(int)
plants["visitor_band"] = pd.cut(plants.wet_visitor_richness,[0,1,4,15,np.inf],labels=False)
plants.to_csv(out/"species.csv")
rows=[]
for col in ["wet_flower","wet_aquatic_hoverfly","bats","amphibians","reptiles","dragonflies_damselflies","spiders"]:
    if col not in plants:continue
    coverage = plants.flower_db if col=="wet_flower" else plants.wet_flower if col=="wet_aquatic_hoverfly" else plants.named_flower
    v = plants[col].astype(float).where(coverage)
    rows.append(dict(metric=col,**compare(plants,v,["site_region","island_band","site_band"],np.random.default_rng(9182))))
for guild, group in chains.groupby("guild"):
    eligible = plants.wet_visitor_richness.gt(0)
    v = plants.accepted.isin(group.plant).astype(float).where(eligible)
    rows.append(dict(metric="wet_degree_adjusted_"+guild,
                     **compare(plants,v,["site_region","island_band","site_band","visitor_band"],np.random.default_rng(9182))))
pd.DataFrame(rows).to_csv(out/"comparisons.csv",index=False)
print("archive",db.sql("select count(*) from globi").fetchone(),"direct",len(direct),"exact",len(same))
