"""Documented-guild and genus commonality of old-only GloBI moth-positive plants.

Pure post-hoc diagnostics on the already frozen 2019/21/23/24/25 same-site
photo cohort: no local actual moth visitation, no proven absent guild.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ID = "wcvp_accepted_plant_name_id"


def bh(ps):
    p = np.asarray(ps, dtype=float)
    if not len(p):
        return np.array([], dtype=float)
    order = np.argsort(p)
    vals = np.minimum.accumulate(
        (p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(vals, 1.0)
    return out


def test_two_props(a: int, n: int, b: int, m: int):
    odds,p = fisher_exact([[a,n-a],[b,m-b]], alternative="two-sided")
    corrected = ((a+.5)*(m-b+.5))/((n-a+.5)*(b+.5))
    return dict(odds_ratio_Haldane=float(corrected), p_nominal=float(p))


def commonality(sp: pd.DataFrame, old: pd.DataFrame) -> dict:
    species = sp.drop_duplicates(ID)[[ID,"wcvp_accepted_name","class"]].copy()
    assert species[ID].is_unique and not species.wcvp_accepted_name.isna().any()
    species["genus"] = species.wcvp_accepted_name.str.split().str[0]
    retained = species[species["class"].isin(["old_only","common"])].copy()
    assert (int(retained["class"].eq("old_only").sum()),
            int(retained["class"].eq("common").sum())) == (713,612)

    tab = pd.crosstab(retained.genus,retained["class"]).reindex(
        columns=["old_only","common"], fill_value=0)
    tab["total"] = tab.old_only + tab.common
    top = tab.sort_values(["old_only","total"],ascending=False).head(25)
    genus_tests = []
    for genus,rec in tab.loc[tab.total.ge(5)].iterrows():
        a,b = int(rec.old_only),int(rec.common)
        record = test_two_props(a,713,b,612)
        genus_tests.append(dict(genus=genus,old_only_species=a,common_species=b,
                                all_old_positive_species=a+b,**record))
    if genus_tests:
        qq=bh([z["p_nominal"] for z in genus_tests])
        for z,q in zip(genus_tests,qq):z["q_BH"]=float(q)
    genus_tests.sort(key=lambda z:(z["q_BH"],z["p_nominal"]))
    enriched=[z for z in genus_tests if z["odds_ratio_Haldane"]>1]
    depleted=[z for z in genus_tests if z["odds_ratio_Haldane"]<1]

    moth = old.loc[old.guild.eq("moth_family_visit") &
                   old.plant_name.isin(retained.wcvp_accepted_name)].copy()
    moth["inat_title"] = moth.study_title.fillna("").str.contains(
        "inaturalist",case=False,regex=False)
    moth["moth_pollinates_claim"] = moth.relation.eq("pollinates")
    origin = moth.groupby("plant_name").agg(
        n_source_rows=("plant_name","size"),
        n_unique_studies=("study_title","nunique"),
        n_distinct_insect_names=("insect_name","nunique"),
        only_inaturalist_titled=("inat_title","all"),
        has_db_pollinates=("moth_pollinates_claim","any"))
    retained = retained.merge(origin,left_on="wcvp_accepted_name",right_index=True,
                              how="left",validate="many_to_one")
    assert retained.n_source_rows.notna().all()
    bt_names = set(old.loc[old.guild.eq("butterfly_family_visit"),"plant_name"])
    retained["has_old_butterfly_claim"] = retained.wcvp_accepted_name.isin(bt_names)
    flags=["only_inaturalist_titled","has_db_pollinates","has_old_butterfly_claim"]
    group_summaries={}
    for name in ["old_only","common"]:
        p=retained[retained["class"].eq(name)]
        group_summaries[name]=dict(
            n_species=len(p),
            n_genera=int(p.genus.nunique()),
            median_relationship_rows_per_species=float(p.n_source_rows.median()),
            median_distinct_study_titles=float(p.n_unique_studies.median()),
            median_documented_insect_names=float(p.n_distinct_insect_names.median()),
            **{f"n_{k}":int(p[k].sum()) for k in flags},
            **{f"fraction_{k}":float(p[k].mean()) for k in flags})
    documentation_tests=[]
    for k in flags:
        a=group_summaries["old_only"][f"n_{k}"]
        b=group_summaries["common"][f"n_{k}"]
        documentation_tests.append(dict(source_documentation_attribute=k,
                                        old_only_positive=a, common_positive=b,
                                        **test_two_props(a,713,b,612)))
    for d,q in zip(documentation_tests,bh([x["p_nominal"] for x in documentation_tests])):
        d["q_BH_3"]=float(q)

    return dict(
        design="713 old-only versus 612 retained old-positive WCVP accepted photographed plant species, old bounded GloBI study titles; NO LOCAL POLLINATION",
        group_summaries=group_summaries,
        n_genera_both_positive_classes=int(retained.genus.nunique()),
        n_genera_shared_between_713_and_612=int(((tab.old_only>0)&(tab.common>0)).sum()),
        top25_old_only_genera=[
            dict(genus=str(g),old_only_species=int(row.old_only),
                 common_species=int(row.common)) for g,row in top.iterrows()],
        genus_tests_over_min5_species=int(len(genus_tests)),
        genus_old_only_enrichment_tests_sorted_q=enriched[:20],
        genus_old_only_depletion_tests_sorted_q=depleted[:20],
        source_documentation_comparison= documentation_tests,
        caveats=["Butterfly documentation is not effective diurnal compensation",
                 "Database pollinates labels are not verified pollen deposition",
                 "No-source record is unknown, not true biological absence",
                 "Multiple genus tests post-hoc; taxa phylogenetically dependent",
                 "Source study title is an approximation of evidence provenance"]
    )
