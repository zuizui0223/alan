"""Attribute the existing frozen fit to photographed species, without taxon-specific model selection."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from alan_gecko_primary_fiveyear import primary
from alan_globi_v011_same_fiveyear import ID

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"


def influence(sites):
    # Raw-response FWL influence, as in PR #43. Keep both light coefficients.
    cols = ["island", "loglight", "light_island", "logpop", "logphotos",
            "logspecies", "bio5", "lograin", "elevation_km"]
    groups = sites.region_year
    n = sites.groupby(["region_year", "land_class"]).key.transform("size")
    weights = n.groupby(groups).transform("min") / n
    x = sites[cols]
    mu = x.mul(weights, axis=0).groupby(groups).transform("sum").div(weights.groupby(groups).transform("sum"),axis=0)
    centered = (x-mu).mul(np.sqrt(weights),axis=0).to_numpy()
    return np.linalg.solve(centered.T @ centered,centered.T)[[1,2],:] * np.sqrt(weights.to_numpy())


def summarize_morphology(cohort, path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != "7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577":
        raise ValueError("Expected the source-direct trait release used by PR 44")
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    axes = ["flower_primary_color","floral_form","tube_depth_class"]
    d = d[d.trait_name.isin(axes) & d.quality.isin(["high","medium"]) & d.resolution_status.eq("resolved")]
    assert not d.duplicated(["accepted_species","trait_name"]).any()
    d = d[d.multistate.eq("False") & d.normalized_value.ne("") & ~d.normalized_value.str.contains("|",regex=False)]
    a = d.pivot(index="accepted_species",columns="trait_name",values="normalized_value").reindex(cohort.index)
    broad_known = a[axes[:2]].notna().all(axis=1)
    broad = broad_known & a.flower_primary_color.eq("white") & a.floral_form.isin(["tubular","salverform","funnel_trumpet"])
    strict_known = a[axes].notna().all(axis=1)
    summary = pd.DataFrame(dict(source_class=cohort.source_class, species=1,broad_known=broad_known,
                                broad_candidate=broad,strict_known=strict_known,
                                strict_candidate=broad & strict_known & a.tube_depth_class.eq("deep")))
    # Only counts by the previously defined GloBI strata are public.
    return summary.groupby("source_class").sum().astype(int)


def run(inputs, signed, out, trait_ledger):
    sites, digest, members, _ = primary(inputs,signed)
    traits = pd.read_csv(RESULTS/"moth_pollen_evidence/shared_544_evidence.csv").set_index("plant")
    hosts = members.rename(columns={"wcvp_accepted_name":"plant"}).copy()
    flags = pd.read_csv(RESULTS/"globi_v011_fiveyear/source_moth_species.csv",dtype={"accepted_id":str})
    hosts["v011"] = hosts[ID].astype(str).isin(flags.loc[flags.v011_wcvp,"accepted_id"])
    assert not hosts.duplicated(["key",ID]).any()
    assert np.array_equal(sites.set_index("key").n_species.sort_index(),hosts.groupby("key").size())
    slopes = influence(sites)
    design = sites[["key","land_class","year","region5","region_year"]].copy()
    design["mainland_pp"] = 100*slopes[0]/sites.n_species
    design["gap_pp"] = 100*slopes[1]/sites.n_species
    obs = hosts.merge(design,on="key",validate="m:1")
    by = obs.groupby("plant")
    estimates = by.agg(accepted_id=(ID,"first"),legacy=("all_moth","first"),v011=("v011","first"),
                       mainland_pp=("mainland_pp","sum"),gap_pp=("gap_pp","sum"))
    estimates["island_pp"] = estimates.mainland_pp + estimates.gap_pp
    assert by[ID].nunique().max()==1
    np.testing.assert_allclose(estimates.loc[estimates.legacy,"gap_pp"].sum(),-1.901720622055596,atol=1e-9)
    np.testing.assert_allclose(estimates.loc[estimates.v011,"gap_pp"].sum(),0.750257336,atol=1e-8)
    ledger = traits.join(estimates,how="left",validate="1:1")
    assert len(ledger)==544 and ledger.legacy.all() and ledger.accepted_id.notna().all()
    ledger["source_class"] = np.where(ledger.v011,"retained_in_v011","old_only")
    measured = obs[obs.plant.isin(ledger.index)]
    count = measured.groupby(["plant","land_class"]).key.nunique().unstack()
    ledger = ledger.join(count.rename(columns={"island":"island_site_years","mainland":"mainland_site_years"}))
    pairs = measured.groupby(["plant","region_year"]).land_class.nunique().eq(2)
    ledger["region_years_with_species_on_both_sides"] = pairs.groupby("plant").sum()
    ledger["years_observed"] = measured.groupby("plant").year.agg(lambda x:"|".join(map(str,sorted(x.unique()))))
    photos = pd.read_csv(RESULTS/"phenology_support/species_544.csv").set_index("plant")
    ledger = ledger.join(photos[["photos","regions","flowering_annotations"]],validate="1:1")

    review = pd.read_csv(RESULTS/"pl_transport_species_review/review_53.csv").set_index("plant")
    promoted = review[review.status.eq("confirmed")]
    ledger["transport_review_status"] = np.where(ledger.moth_transport_documented | ledger.index.isin(promoted.index),
                                                  "documented","unconfirmed")
    original = pd.read_csv(RESULTS/"moth_pollen_evidence/pollen_evidence.csv")
    ledger["transport_evidence"] = original.groupby("plant").measurement.agg(lambda x:" | ".join(sorted(set(x))))
    ledger["transport_source"] = original.groupby("plant").source.agg(lambda x:" | ".join(sorted(set(x))))
    ledger.loc[promoted.index,"transport_evidence"] = promoted.reason
    ledger.loc[promoted.index,"transport_source"] = promoted.source
    night = pd.read_csv(RESULTS/"night_exclusion/species.csv").set_index("plant")
    np.testing.assert_allclose(ledger.loc[night.index,"gap_pp"],night.gap,atol=1e-10)
    np.testing.assert_array_equal(ledger.loc[night.index,"island_site_years"],night.island_positive_site_years)
    ledger = ledger.join(night[["loss","studies"]].rename(columns={"loss":"day_only_relative_loss","studies":"night_exclusion_studies"}))
    pl = ledger.dropna(subset="glopl_median_log_response_ratio")
    assert len(pl)==64 and pl.transport_review_status.eq("documented").sum()==17
    np.testing.assert_allclose(np.polyfit(pl.glopl_median_log_response_ratio,pl.gap_pp,1)[0],-0.0193800457560167,atol=1e-12)
    ledger["direct_flower_form_status"] = "species_label_not_in_public_trait_results"
    ledger["difference_if_relabelled_v011_pp"] = np.where(ledger.v011,0,-ledger.gap_pp)
    out.mkdir(parents=True,exist_ok=True)
    ledger.to_csv(out/"shared_544.csv")
    summarize_morphology(ledger,trait_ledger).to_csv(out/"morphology_counts.csv")
    print(digest,ledger.groupby("source_class").gap_pp.agg(["size","sum"]).to_dict())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("inputs","signed","output","trait_ledger"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    run(a.inputs,a.signed,a.output,a.trait_ledger)
