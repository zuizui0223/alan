"""Build private compact GBIF inputs from the published Info Flora Darwin Core archive."""
import sys
from pathlib import Path
from zipfile import ZipFile
import duckdb
import pyarrow as pa
import pyarrow.csv as csv
import pyarrow.parquet as pq

columns="""countryCode basisOfRecord occurrenceStatus coordinateUncertaintyInMeters decimalLatitude decimalLongitude eventDate recordedBy genus specificEpithet phylum taxonRank lifeStage bibliographicCitation references associatedReferences associatedMedia""".split()
def convert(archive,cache):
    cache.mkdir(parents=True,exist_ok=True)
    database=duckdb.connect(config={"threads":3})
    out=cache/"eligible.parquet"
    writer=None
    n=0
    try:
        with ZipFile(archive) as z, z.open("occurrence.txt") as src:
            for batch in csv.open_csv(src,read_options=csv.ReadOptions(block_size=32*1024*1024),
                                    parse_options=csv.ParseOptions(delimiter="\t",quote_char=False),
                                    convert_options=csv.ConvertOptions(include_columns=columns,column_types={c:pa.string() for c in columns})):
                n+=len(batch)
                database.register("records",pa.Table.from_batches([batch]))
                rows=database.sql("""
                  select try_cast(eventDate as date) as day, decimalLatitude::double lat, decimalLongitude::double lon,
                    trim(recordedBy) as recorder, genus||' '||specificEpithet as species, lifeStage
                  from records where countryCode='CH' and lower(basisOfRecord)='humanobservation'
                    and lower(occurrenceStatus)='present' and phylum='Tracheophyta'
                    and try_cast(coordinateUncertaintyInMeters as double) between 0 and 3536
                    and try_cast(decimalLatitude as double) between 45 and 48.5
                    and try_cast(decimalLongitude as double) between 5 and 11
                    and regexp_full_match(eventDate,'[0-9]{4}-[0-9]{2}-[0-9]{2}')
                    and try_cast(eventDate as date) between date '1992-01-01' and date '2024-12-31'
                    and regexp_full_match(genus,'[A-Z][a-z-]+') and regexp_full_match(specificEpithet,'[a-z][a-z-]+')
                    and lower(taxonRank) in ('species','subspecies','variety','varietas','form','forma')
                    and not contains(lower(concat_ws(' ',bibliographicCitation,"references",associatedReferences,associatedMedia)), 'inaturalist')
                """).arrow().read_all()
                if writer is None:writer=pq.ParquetWriter(out.with_suffix(".partial"),rows.schema,compression="zstd")
                writer.write_table(rows)
    finally:
        if writer:writer.close()
    assert n==16302425,n
    out.with_suffix(".partial").replace(out)
    print(n,flush=True)

if __name__=="__main__":convert(*map(Path,sys.argv[1:3]))

