"""Extract named moth prey of the two GBIF house geckos from the frozen release."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as csv

from alan_globi_predator_pollinator_multilayer_v011 import (
    FIELDS, SOURCE_ROWS, SOURCE_SHA, annotate, sha256
)

GECKOS = ("Hemidactylus frenatus", "Hemidactylus mabouia")


def extract(source, destination):
    assert sha256(source) == SOURCE_SHA
    with gzip.open(source, "rt", encoding="utf-8") as handle:
        import csv as std_csv
        columns = set(next(std_csv.reader(handle)))
    included = sorted(columns & FIELDS)
    names = [side + field for side in ("source", "target") for field in ("TaxonSpeciesName", "TaxonName")
             if side + field in included]
    records = []
    n = 0
    with pa.input_stream(str(source), compression="gzip") as stream:
        batches = csv.open_csv(stream, read_options=csv.ReadOptions(block_size=8 << 20),
            parse_options=csv.ParseOptions(newlines_in_values=True),
            convert_options=csv.ConvertOptions(include_columns=included,
                                                column_types={k: pa.string() for k in included}))
        for batch in batches:
            n += batch.num_rows
            mask = pc.is_in(batch.column(names[0]), value_set=pa.array(GECKOS))
            for col in names[1:]:
                mask = pc.or_(mask, pc.is_in(batch.column(col), value_set=pa.array(GECKOS)))
            found = batch.filter(mask)
            if not found.num_rows:
                continue
            for group, actor, role, prey, other_group, stage, scopes in annotate(found.to_pandas().astype("string")):
                if group == "geckos" and actor in GECKOS and role == "moth_predation":
                    records.append((actor, prey, stage, "non_inaturalist_trophich" in scopes))
    assert n == SOURCE_ROWS and set(r[0] for r in records) == set(GECKOS)
    destination.mkdir(exist_ok=True, parents=True)
    rows = pd.DataFrame(records,columns=["gecko","moth","source_stage","non_inat_trophich"])
    out = destination / "species_prey.csv"
    rows.value_counts(sort=False).rename("relationship_rows").reset_index().to_csv(out,index=False)
    (destination / "receipt.json").write_text(json.dumps(dict(scanned=n, source_sha256=SOURCE_SHA,
        rows=len(records), prey_table_sha256=hashlib.sha256(out.read_bytes()).hexdigest()),indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source",type=Path)
    p.add_argument("destination",type=Path)
    args=p.parse_args()
    extract(args.source,args.destination)
