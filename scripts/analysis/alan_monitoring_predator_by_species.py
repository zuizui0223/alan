"""Run the existing monitoring/predator overlap separately for each four-arm plant."""
import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from alan_repeat_flower_contrasts import four, GROUP, read_source
from alan_repeat_flower_predators import run


def summarize(work, output):
    folder = work / "phenology_run"
    sources = [("NPN", "npn/individual_visits.parquet", "npn_exposure.parquet"),
               ("Sweden", "sweden_all_visits.parquet", "sweden_all_exposure.parquet"),
               ("ClimateWatch", "australia_standard_visits.parquet", "australia_exposure.parquet")]
    observations = pd.concat([read_source(name, folder / visits, folder / light)
                              for name, visits, light in sources], ignore_index=True)
    candidates = four(observations).index.to_frame(index=False)
    assert len(candidates) == 5
    subset = observations.merge(candidates, on=GROUP, validate="m:1")
    rows = []
    with TemporaryDirectory() as tmp:
        path = Path(tmp)
        for plant, group in subset.groupby("plant"):
            run(group, work / "predator", work / "v011/coast", folder / "pred_days", path)
            summary = json.loads((path / "predator_support.json").read_text())
            pairs = pd.read_csv(path / "predator_overlap.csv")
            rows.append(dict(plant=plant, n_status=len(group), n_open=int(group.status.sum()),
                             **summary, predator_labels="|".join(sorted(pairs.predator.unique()))))
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    summarize(args.work, args.output)
