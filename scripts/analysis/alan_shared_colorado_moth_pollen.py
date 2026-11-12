"""Reconcile species-resolved moth pollen with the frozen 544, using a pinned Dryad version."""
import argparse
import hashlib
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
VERSION = "https://datadryad.org/api/v2/versions/329747/download"
WORKBOOK = "moths_cleaner2.xlsx"
SHA256 = "b9c99f6a6d92b926c632290256290688a11738892a8366e4832cdca54cc8ed73"


def summarize(archive, output):
    if not archive.exists():
        response = requests.get(VERSION, timeout=90)
        response.raise_for_status()
        archive.parent.mkdir(parents=True, exist_ok=True)
        archive.write_bytes(response.content)
    with ZipFile(archive) as files:
        workbook = files.read(WORKBOOK)
    if hashlib.sha256(workbook).hexdigest() != SHA256:
        raise ValueError("Dryad workbook differs from archived file 3652358")
    data = pd.read_excel(BytesIO(workbook), sheet_name="moths_cleaner2")
    assert len(data) == data["Individual #"].nunique() == 276 and data["Recap?"].eq("N").all()
    cohort = pd.read_csv(ROOT / "results/species_interpretation/shared_544.csv").set_index("plant")
    assert len(cohort) == 544 and cohort.transport_review_status.eq("documented").sum() == 49
    region = pd.read_csv(ROOT / "results/bilateral_visitors/regional_frame.csv")
    assert len(region) == 195 and region.region5.nunique() == 78 and not region.region5.eq("14_25").any()
    pollen_columns = data.columns[data.columns.get_loc("Significant pollen load") + 1:]
    matched = pollen_columns.intersection(cohort.index)
    assert len(matched) == 3
    assert cohort.loc[matched, "transport_review_status"].eq("unconfirmed").all()
    assert cohort.loc[matched, "glopl_cases"].eq(0).all()
    seasons = []
    edges = []
    for plant in matched:
        for season, group in data.groupby("Season"):
            measured = group[plant]
            carried = group[measured.ge(5)]
            assert measured.notna().all() and set(measured[measured.gt(0)]) == set(carried[plant])
            seasons.append(dict(plant=plant, season=season, moths_examined=len(group),
                                carriers=len(carried), grains=int(carried[plant].sum()),
                                source_review=cohort.at[plant, "transport_review_status"],
                                source_sha256=SHA256))
            for label, insects in carried.groupby("Species", dropna=False):
                edges.append(dict(plant=plant, season=season, original_moth_label=label,
                                  carriers=len(insects), grains=int(insects[plant].sum())))
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(seasons).to_csv(output / "colorado_pollen.csv", index=False)
    pd.DataFrame(edges).to_csv(output / "colorado_pollen_carriers.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    summarize(args.archive, args.output)
