from alan_extended11 import coords,get,batch
from alan_decadal_light import C
def phen(job):
    st,start=job
    r=coords.loc[st]
    path=C/f"green_{st}_{start}.json"
    if path.exists():return
    try:
        d=get(path,"https://modis.ornl.gov/rst/api/v1/MCD12Q2/subset",dict(latitude=r.latitude,longitude=r.longitude,
            startDate=f"A{start}001",endDate=f"A{min(start+7,2024)}001",kmAboveBelow=0,kmLeftRight=0))
        assert d["nrows"]==d["ncols"]==1
    except Exception as e: print(st,start,str(e)[:50],flush=True)
if __name__=="__main__":batch(phen,[(s,y) for s in coords.index for y in (2001,2009,2017)],12)
