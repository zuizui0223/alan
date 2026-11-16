"""Read JMA on-premises first flowering, using the station position on the event date."""
from pathlib import Path
from zipfile import ZipFile
import io
import sys
import pandas as pd
import requests

# Single-species normal targets; omit aggregated cherry, dandelion and violet protocols.
TAXA={"052":"Trifolium repens","090":"Syringa vulgaris","088":"Kerria japonica",
      "085":"Prunus persica","030":"Prunus armeniaca","037":"Platycodon grandiflorus",
      "010":"Lagerstroemia indica","007":"Wisteria floribunda","073":"Lycoris radiata",
      "063":"Lilium longiflorum","011":"Miscanthus sinensis","035":"Diospyros kaki",
      "042":"Castanea crenata"}


def read_jma(root):
    root.mkdir(exist_ok=True)
    base="https://www.data.jma.go.jp/sakura/data/"
    from bs4 import BeautifulSoup
    from urllib.parse import urljoin
    page=urljoin(base,"download_ruinenchi.html")
    if not (root/"all.zip").exists():
        r=requests.get(page,timeout=30);r.raise_for_status()
        archive=next(a["href"] for a in BeautifulSoup(r.content,"html.parser").select("a[href]") if ".zip" in a["href"])
        r=requests.get(urljoin(page,archive),timeout=45);r.raise_for_status()
        (root/"all.zip").write_bytes(r.content)
    if not (root/"smaster.index.zip").exists():
        r=requests.get("https://www.data.jma.go.jp/stats/data/mdrr/chiten/meta/smaster.index.zip",timeout=45)
        r.raise_for_status();(root/"smaster.index.zip").write_bytes(r.content)
    with ZipFile(root/"smaster.index.zip") as z:
        lines=z.read("smaster.index").splitlines()
    history=[]
    for line in lines:
        lat=line[36:42];lon=line[42:49]
        if not (lat.isdigit() and lon.isdigit()):continue
        history.append((int(line[:3]),int(line[64:72]),int(line[72:80].strip() or b"99999999"),
                        int(lat[:2])+int(lat[2:])/6000,int(lon[:3])+int(lon[3:])/6000))
    history=pd.DataFrame(history,columns=["station","first","last","latitude","longitude"])
    rows=[]
    with ZipFile(root/"all.zip") as z:
        for code,species in TAXA.items():
            f=pd.read_csv(io.BytesIO(z.read(f"ruinenchi_all/{code}.csv")),encoding="cp932",skiprows=1)
            for year in range(1992,2015):
                ix=f.columns.get_loc(str(year))
                # 8 = designated species on station grounds; 6/7 are alternative taxa.
                q=f.iloc[:,[0,ix,ix+1]].copy();q.columns=["station","mmdd","flag"]
                q=q[q.flag.eq(8)&q.mmdd.gt(0)].copy()
                q["date"]=year*10000+q.mmdd
                q=q.merge(history,on="station").query("first<=date<=last")
                q["flower"]=pd.to_datetime(q.date.astype(int).astype(str),format="%Y%m%d").dt.dayofyear
                rows.append(q.assign(year=year,species=species,source="JMA")[["station","year","latitude","longitude","flower","species","source"]])
    d=pd.concat(rows).drop_duplicates()
    assert not d.duplicated(["station","species","year"]).any()
    return d


if __name__=="__main__":
    d=read_jma(Path(sys.argv[1]));print(d.groupby("species").size())
