"""Estimate species-level slope uncertainty on both sides of the frozen ALAN model."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, t

from alan_gecko_primary_fiveyear import primary

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
COVARS = ["island", "loglight", "light_island", "logpop", "logphotos", "logspecies",
          "bio5", "lograin", "elevation_km"]


def centered(values, groups, weights):
    sums = np.zeros((groups.max() + 1, values.shape[1]))
    np.add.at(sums, groups, values * weights[:, None])
    means = sums / np.bincount(groups, weights=weights)[:, None]
    return (values - means[groups]) * np.sqrt(weights[:, None])


def run(inputs, signed, out):
    sites, _, members, _ = primary(inputs, signed)
    traits = pd.read_csv(RESULTS / "species_interpretation/shared_544.csv").set_index("plant")
    assert len(sites) == 10698 and sites.key.is_unique and len(traits) == 544
    sites = sites.reset_index(drop=True)
    cells = pd.factorize(sites.region_year)[0]
    regions = pd.factorize(sites.region5)[0]
    nregions = len(np.unique(regions))
    assert nregions == 78 and len(np.unique(cells)) == 195
    counts = sites.groupby(["region_year", "land_class"]).key.transform("size")
    weights = (counts.groupby(sites.region_year).transform("min") / counts).to_numpy()
    x = centered(sites[COVARS].to_numpy(), cells, weights)
    photos = members[members.wcvp_accepted_name.isin(traits.index)]
    row = pd.Index(sites.key).get_indexer(photos.key)
    col = traits.index.get_indexer(photos.wcvp_accepted_name)
    assert min(row) >= 0 and not photos.duplicated(["key", "wcvp_accepted_name"]).any()
    y = np.zeros((len(sites), len(traits)))
    y[row, col] = 100 / sites.n_species.to_numpy()[row]
    y = centered(y, cells, weights)
    bread = np.linalg.inv(x.T @ x)
    beta = bread @ x.T @ y
    residual = y - x @ beta
    mult = nregions / (nregions - 1) * (len(x) - 1) / (len(x) - len(COVARS))
    for side, contrast in [("mainland", bread[1]), ("island", bread[1] + bread[2]), ("gap", bread[2])]:
        point = contrast @ x.T @ y
        np.testing.assert_allclose(point, traits[side + "_pp"], atol=1e-10)
        cluster_scores = np.zeros((nregions, len(traits)))
        np.add.at(cluster_scores, regions, (x @ contrast)[:, None] * residual)
        se = np.sqrt(mult * (cluster_scores ** 2).sum(axis=0))
        assert np.all(se > 0)
        half = t.ppf(.975, nregions - 1) * se
        traits[side + "_lo"] = point - half
        traits[side + "_hi"] = point + half
        traits[side + "_p"] = 2 * t.sf(abs(point / se), nregions - 1)
    tests = traits[["island_p", "mainland_p"]].to_numpy()
    adjusted = false_discovery_control(tests.ravel()).reshape(tests.shape)
    traits["island_q"], traits["mainland_q"] = adjusted.T
    traits["island_rank"] = traits.island_pp.rank().astype(int)
    traits["mainland_rank"] = traits.mainland_pp.rank().astype(int)
    extra = pd.read_csv(RESULTS / "mechanism_exemplars/transport_updates.csv").set_index("plant")
    traits.loc[extra.index, "transport_review_status"] = "documented"
    traits.loc[extra.index, "transport_evidence"] = extra.evidence
    links = pd.read_csv(RESULTS / "dual_role_80/island_cohort_reviewed_links.csv")
    pred = links[links.assessment.eq("predation_recorded")]
    traits["reviewed_dual_role_predator"] = traits.index.map(pred.groupby("plant").predator.first())
    traits["reviewed_dual_role_adult_prey"] = traits.index.map(pred.groupby("plant").moth.first())
    fields = ["island_pp", "island_lo", "island_hi", "island_q", "island_rank", "mainland_pp",
              "mainland_lo", "mainland_hi", "mainland_q", "mainland_rank", "gap_pp", "gap_lo", "gap_hi",
              "island_site_years", "mainland_site_years", "region_years_with_species_on_both_sides",
              "v011", "transport_review_status", "transport_evidence", "day_only_relative_loss",
              "glopl_cases", "glopl_median_log_response_ratio", "sc_source_class",
              "reviewed_dual_role_predator", "reviewed_dual_role_adult_prey"]
    eligible = traits.island_pp.lt(0) | traits.mainland_pp.lt(0)
    traits.loc[eligible, fields].to_csv(out)
    for side in ["island", "mainland"]:
        print(side, "negative", sum(traits[side+"_pp"] < 0), "nominal negative", sum(traits[side+"_hi"] < 0),
              "BH negative", sum((traits[side+"_pp"] < 0) & (traits[side+"_q"] < .05)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs", "signed", "output"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    run(args.inputs, args.signed, args.output)
