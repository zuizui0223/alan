"""Project public GBIF snapshot columns into a local, resumable worldwide occurrence database."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import time
import json
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.fs as fs
import pyarrow.parquet as pq

SNAPSHOT="gbif-open-data-eu-central-1/occurrence/2026-11-01/occurrence.parquet"
COLUMNS=["specieskey","species","datasetkey","basisofrecord","occurrencestatus","decimallatitude","decimallongitude",
         "coordinateuncertaintyinmeters","year","month","day","issue"]
ACCEPT=["HUMAN_OBSERVATION","OBSERVATION","PRESERVED_SPECIMEN"]
def scan(folder,workers,start,end,annual=False):
    folder.mkdir(exist_ok=True,parents=True)
    s3=fs.S3FileSystem(anonymous=True,region="eu-central-1",request_timeout=120,connect_timeout=20)
    objects=sorted(x.base_name for x in s3.get_file_info(fs.FileSelector(SNAPSHOT)) if x.base_name.isdigit())
    assert len(objects)==10346
    def process(key):
        done=folder/(key+".json")
        if done.exists():return
        for attempt in range(5):
            try:
                with s3.open_input_file(SNAPSHOT+"/"+key) as fp:
                    source=pq.ParquetFile(fp,pre_buffer=False)
                    n=source.metadata.num_rows
                    plant=source.read(columns=["phylum"],use_threads=False)["phylum"]
                    mask=pc.fill_null(pc.equal(plant,"Tracheophyta"),False)
                    count=pc.sum(mask).as_py() or 0
                    kept=0
                    if count:
                        tab=source.read(columns=COLUMNS,use_threads=False).filter(mask)
                        y=tab["year"]
                        mask=pc.or_(pc.and_(pc.greater_equal(y,1992),pc.less_equal(y,1999)),
                                    pc.and_(pc.greater_equal(y,2020),pc.less_equal(y,2024)))
                        if annual:
                            mask=pc.and_(pc.greater_equal(y,1992),pc.less_equal(y,2024))
                        def add(v):
                            nonlocal mask
                            mask=pc.and_(mask,v)
                        add(pc.is_valid(tab["specieskey"]))
                        add(pc.is_valid(tab["species"]))
                        add(pc.equal(tab["occurrencestatus"],"PRESENT"))
                        add(pc.is_in(tab["basisofrecord"],value_set=pa.array(ACCEPT)))
                        add(pc.greater_equal(tab["coordinateuncertaintyinmeters"],0))
                        add(pc.less_equal(tab["coordinateuncertaintyinmeters"],1000))
                        for col,a,b in [("decimallongitude",-180,180),("decimallatitude",-90,90),("day",1,31),("month",1,12)]:
                            add(pc.greater_equal(tab[col],a));add(pc.less_equal(tab[col],b))
                        tab=tab.filter(mask)
                        issues=tab["issue"]
                        if len(tab):
                            invalid={"ZERO_COORDINATE","COORDINATE_OUT_OF_RANGE","COORDINATE_INVALID","COUNTRY_COORDINATE_MISMATCH"}
                            if pa.types.is_list(issues.type) or pa.types.is_large_list(issues.type):
                                vals=pc.list_flatten(issues)
                                if pa.types.is_struct(vals.type):vals=pc.struct_field(vals,"array_element")
                                badrows=pc.filter(pc.list_parent_indices(issues),pc.is_in(vals,value_set=pa.array(sorted(invalid)))).to_numpy()
                                valid=np.ones(len(tab),bool);valid[badrows]=False
                                tab=tab.filter(pa.array(valid))
                            else:
                                tab=tab.filter(pc.invert(pc.match_substring_regex(pc.fill_null(issues,""),"|".join(sorted(invalid)))))
                        tab=tab.drop(["issue","occurrencestatus","coordinateuncertaintyinmeters"])
                        kept=len(tab)
                        if kept:
                            temp=folder/(key+".tmp")
                            pq.write_table(tab,temp,compression="zstd");temp.replace(folder/(key+".parquet"))
                    done.with_suffix(".part").write_text(json.dumps([n,count,kept]))
                    done.with_suffix(".part").replace(done)
                return
            except Exception:
                if attempt==4:raise
                time.sleep(2**attempt)
    with ThreadPoolExecutor(workers) as pool:
        for _ in pool.map(process,objects[start:end]):pass

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("folder",type=Path);parser.add_argument("--workers",type=int,default=12)
    parser.add_argument("--start",type=int,default=0);parser.add_argument("--end",type=int,default=None)
    parser.add_argument("--annual", action="store_true")
    args=parser.parse_args();scan(args.folder,args.workers,args.start,args.end,args.annual)

