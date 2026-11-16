"""Contrast trait composition against other screened plants within region and support band."""
from pathlib import Path
from zipfile import ZipFile
import sys
import numpy as np
import pandas as pd


def run(wcvp,ledger,out):
    data=pd.read_csv(out/"all_regions.csv.gz")
    names=set(data.species)
    cols=["taxon_name","plant_name_id","accepted_plant_name_id","taxon_status","taxon_rank","family","lifeform_description"]
    with ZipFile(wcvp) as z:
        records=pd.read_csv(z.open("wcvp_names.csv"),sep="|",usecols=cols,dtype=str,low_memory=False)
    accepted=records[records.taxon_status.eq("Accepted")&records.taxon_rank.str.lower().eq("species")].copy()
    aliases=records[records.taxon_name.isin(names)&records.taxon_status.isin(["Accepted","Synonym"])].copy()
    aliases["id"]=aliases.accepted_plant_name_id.where(~aliases.taxon_status.eq("Accepted"),aliases.plant_name_id)
    aliases=aliases[aliases.groupby("taxon_name").id.transform("nunique").eq(1)].drop_duplicates("taxon_name")
    traits=aliases.merge(accepted[["plant_name_id","taxon_name","family","lifeform_description"]],left_on="id",right_on="plant_name_id",suffixes=("","_accepted"))
    traits=traits.rename(columns={"taxon_name":"species","taxon_name_accepted":"accepted_name","family_accepted":"wcvp_family","lifeform_description_accepted":"lifeform"})
    e=pd.read_csv(ledger).rename(columns={"plant":"accepted_name"})
    traits=traits[["species","accepted_name","wcvp_family","lifeform"]].merge(e[["accepted_name","transport_review_status","sc_source_class","glopl_median_log_response_ratio","day_only_relative_loss"]],how="left",validate="m:1")
    desc=traits.lifeform
    for key,pat in {"annual":"\\bannual\\b","woody":"tree|shrub|liana","geophyte":"geophyte","perennial_word":"perennial"}.items():
        traits[key]=desc.str.contains(pat,case=False).where(desc.notna()).astype(float)
    traits["moth_transport_documented"] =traits.transport_review_status.eq("documented").where(traits.transport_review_status.isin(["documented","unconfirmed"])).astype(float)
    traits["SI_among_SC_SI"]=traits.sc_source_class.eq("SI").where(traits.sc_source_class.isin(["SI","SC"])).astype(float)
    traits["pollen_limitation"]=traits.glopl_median_log_response_ratio
    traits["night_exclusion_loss"]=traits.day_only_relative_loss
    base=data[data.blocks.ge(3)].merge(traits,how="left",validate="m:1")
    base["band"]=pd.cut(base.blocks,[2,5,10,25,100000],labels=False)
    base[["sp","species","accepted_name","wcvp_family","lifeform","sc_source_class","transport_review_status","pollen_limitation","night_exclusion_loss"]].drop_duplicates().to_csv(out/"species_traits.csv",index=False)
    keys=["annual","woody","geophyte","perennial_word","moth_transport_documented","SI_among_SC_SI","pollen_limitation","night_exclusion_loss"]
    families=base.drop_duplicates("sp").wcvp_family.value_counts()
    for fam in families[families>=30].index:
        col="family_"+fam
        base[col]=base.wcvp_family.eq(fam).where(base.wcvp_family.notna()).astype(float);keys.append(col)
    rows=[]
    for key in keys:
        known=base.dropna(subset=[key])
        for land in ["both","island","mainland"]:
            g=known if land=="both" else known[known.land_class.eq(land)]
            for sign in ["neg","pos"]:
                a=g[g.label.eq(sign)].copy()
                controls=g[~g.label.eq(sign)].groupby(["region","band"])[key].mean().rename("reference")
                a=a.join(controls,on=["region","band"]).dropna(subset=["reference"])
                if a.empty:continue
                diff=a[key]-a.reference
                row=dict(trait=key,land=land,label=sign,known_species=g.sp.nunique(),candidate_species=a.sp.nunique(),candidate_pairs=len(a),regions=a.region.nunique(),
                         observed=a[key].mean(),reference=a.reference.mean(),difference=diff.mean())
                region=diff.groupby(a.region).agg(["sum","size"])
                loo=(diff.sum()-region["sum"])/(len(diff)-region["size"]).replace(0,np.nan)
                row.update(leave_region_min=loo.min(),leave_region_max=loo.max())
                rows.append(row)
    pd.DataFrame(rows).to_csv(out/"trait_contrasts.csv",index=False)

if __name__=="__main__":run(*map(Path,sys.argv[1:4]))
