"""Great Britain atlas baseline-recorded hectads and subsequent re-recording."""
import sys
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

cache, light, out = map(Path, sys.argv[1:4])
out.mkdir(exist_ok=True, parents=True)
atlas = cache / "PA2020_10by10kmSquareData"
names = pd.read_csv(atlas / "nameLookup.csv")
names = names[names["rank"].eq("species") & names.hybrid.isna() & (names.qualifier.isna() | names.qualifier.eq("s.s."))].set_index("ddbidLookup")
pieces = []
with ZipFile(atlas / "hectadDistributions.zip") as archive:
    for file in archive.namelist():
        key = Path(file).stem
        if key not in names.index:
            continue
        rows = pd.read_csv(archive.open(file), usecols=["hectad","1987 - 1999","2010 - 2019"])
        rows.columns = ["cell","early","late"]
        # OSGB National Grid only. Irish and Channel Island systems differ.
        rows = rows[rows.cell.str.fullmatch(r"[HNOST][A-HJ-Z]\d\d") & (rows.early.eq(1) | rows.late.eq(1))].copy()
        rows["id"] = key
        pieces.append(rows)
obs = pd.concat(pieces, ignore_index=True)
assert not obs.duplicated(["cell","id"]).any()
cells = obs.groupby("cell")[["early","late"]].sum()
cells = cells[(cells.early >= 50) & (cells.late >= 50)].copy()
alphabet = "ABCDEFGHJKLMNOPQRSTUVWXYZ"
def origin(cell):
    a,b = [alphabet.index(letter) for letter in cell[:2]]
    return (((a-2)%5*5+b%5)*100000+int(cell[2])*10000,
            (19-a//5*5-b//5)*100000+int(cell[3])*10000)
xy = np.array([origin(c) for c in cells.index])
offsets = np.array([(x,y) for x in range(1000,10000,2000) for y in range(1000,10000,2000)])
samples = xy[:,None,:] + offsets
lon,lat = Transformer.from_crs(27700,4326,always_xy=True).transform(samples[:,:,0].ravel(),samples[:,:,1].ravel())
annual = {}
for year in [*range(1992,2000),*range(2010,2020)]:
    sensor = "calDMSP" if year <= 2013 else "simVIIRS"
    with rasterio.open(light / f"Harmonized_DN_NTL_{year}_{sensor}.tif") as src:
        rr,cc = rasterio.transform.rowcol(src.transform,lon,lat)
        rr,cc = np.asarray(rr),np.asarray(cc)
        win = rasterio.windows.Window.from_slices((rr.min(),rr.max()+1),(cc.min(),cc.max()+1))
        vals = src.read(1,window=win)[rr-rr.min(),cc-cc.min()]
        assert vals.min() >= 0 and vals.max() <= 63
        annual[year] = vals.reshape(-1,25).mean(axis=1)
cells["light_start"] = np.mean([annual[y] for y in range(1992,2000)],axis=0)
cells["delta"] = np.mean([annual[y] for y in range(2010,2020)],axis=0) - cells.light_start
panel = obs[obs.early.eq(1)].join(cells, on="cell",rsuffix="_n",how="inner")
result = []
for key, d in panel.groupby("id"):
    y = d.late.to_numpy(float)
    blocks = d.cell.str[:2]
    nblocks = blocks.nunique()
    row = dict(id=key,species=names.loc[key,"taxonName"],cells=len(d),blocks=nblocks,rerecorded=int(y.sum()))
    result.append(row)
    if len(d)<50 or nblocks<10 or min(sum(y),sum(1-y))<10:
        continue
    # Leave the focal species out of both period-specific richness proxies.
    X = pd.DataFrame(dict(change=d.delta,initial=d.light_start,
                         early_effort=np.log(d.early_n-1),late_effort=np.log(d.late_n-y),response=y))
    within = X-X.groupby(blocks).transform("mean")
    x = within.iloc[:,:4].to_numpy()
    if np.linalg.matrix_rank(x)<4:
        continue
    inv = np.linalg.inv(x.T@x)
    beta = inv @ (x.T@within.response)
    error = within.response.to_numpy()-x@beta
    scores = pd.DataFrame(x*error[:,None],index=d.index).groupby(blocks).sum().to_numpy()
    variance = inv@(scores.T@scores)@inv * nblocks/(nblocks-1)*(len(d)-1)/(len(d)-nblocks-4)
    se = np.sqrt(variance[0,0])
    row["delta_dn_iqr"] = d.delta.quantile(.75)-d.delta.quantile(.25)
    row["pp_per_observed_iqr"] = 100*beta[0]*row["delta_dn_iqr"]
    row.update(pp_per_10dn=1000*beta[0],se_pp=1000*se,p=2*t.sf(abs(beta[0]/se),nblocks-1))
res = names[["taxonName"]].rename(columns={"taxonName":"species"}).join(pd.DataFrame(result).set_index("id").drop(columns="species"))
ok = res.p.notna()
res.loc[ok,"q"] = multipletests(res.loc[ok,"p"],method="fdr_bh")[1]
res.to_csv(out/"atlas_light.csv",index_label="id")
cells.to_csv(out/"grid_support.csv")
print("species, cells, fits",len(res),len(cells),ok.sum())

