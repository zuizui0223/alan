"""Reduce the original signed MOIS register to private observation counts by cell, year and taxon."""
import concurrent.futures
import hashlib
import json
import sys
import tarfile
import tempfile
from pathlib import Path
import duckdb

source, cache = map(Path, sys.argv[1:3])
cache.mkdir(exist_ok=True, parents=True)
receipt = source / "mois_all_original_photo_target_tiers_receipt.json"
assert hashlib.sha256(receipt.read_bytes()).hexdigest() == "4ed5139e510ba8442d4cdb03f2b769d94fac942832b1d544d7a7925bc431255b"
registry = json.loads(receipt.read_text())
def extract(item):
    filename, signature = item
    archive = source / filename
    out = cache / (archive.stem + ".parquet")
    if out.exists():
        return
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == signature, filename
    with tempfile.TemporaryDirectory(dir=cache) as temporary:
        tmp = Path(temporary)
        with tarfile.open(archive) as tar:
            for member in tar:
                if member.isfile() and member.name.endswith(".csv.gz"):
                    (tmp / Path(member.name).name).write_bytes(tar.extractfile(member).read())
        con = duckdb.connect(config={"threads":2,"memory_limit":"1GB"})
        files = (tmp / "*.csv.gz").as_posix()
        # UUIDs are assigned to disjoint observation-hash tiers and source taxon parts.
        sql = f"""
        COPY (
        SELECT try_cast(left(observed_on,4) AS SMALLINT) AS year,
          floor((90 + try_cast(latitude AS DOUBLE))*20)::INTEGER AS y,
          floor((180 + try_cast(longitude AS DOUBLE))*20)::INTEGER AS x,
          try_cast(inat_taxon_id AS BIGINT) AS taxon,
          count(DISTINCT observation_uuid) AS n, count(*) AS photos
        FROM read_csv('{files}', all_varchar=true, header=true)
        WHERE valid_observed_date='true' AND valid_geopoint='true'
          AND try_cast(left(observed_on,4) AS INTEGER) BETWEEN 1992 AND 2024
          AND (try_cast(positional_accuracy AS DOUBLE) IS NULL OR try_cast(positional_accuracy AS DOUBLE)<=2500)
          AND try_cast(latitude AS DOUBLE) BETWEEN -90 AND 90
          AND try_cast(longitude AS DOUBLE) BETWEEN -180 AND 180
        GROUP BY ALL) TO '{(tmp/"out.parquet").as_posix()}' (FORMAT PARQUET)"""
        con.execute(sql)
        con.close()
        (tmp/"out.parquet").replace(out)
if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        list(pool.map(extract,registry["published_tier_archive_SHA256"].items()))
