"""JRC v1.5 occurrence: 1984-2024 BACKGROUND, never 2016 annual water.

Download one small native-grid window using HTTP ranges; cache outside repo.
Nearest-water distances refer to detected pixel centres, not surveyed shoreline.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
import requests
import rasterio
from rasterio.windows import from_bounds, Window
from pyproj import Transformer
from scipy.spatial import cKDTree

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
CACHE=ROOT.parent/'gbif_annual/alan_experiment_sources/giavi_landscape'
CACHE.mkdir(parents=True,exist_ok=True)
URL='https://s3.waw4-1.cloudferro.com/swift/v1/global-surface-water/download2024/Aggregated/VER1-5/occurrence/occurrence_0E_50N_v1_5_2024.tif'
clip=CACHE/'occurrence_v15_1984_2024_swiss_subset.tif'
manifest=CACHE/'occurrence_subset_manifest.json'
sites=pd.read_csv(ROOT/'database/alan_experiments/site_context/primary_coordinates.csv').query("study_id=='Giavi2021'")
assert len(sites)==12
# Deliberately generous margin: >5 km in every direction from all site centres.
bbox=[float(sites.lon.min()-.10),float(sites.lat.min()-.08),float(sites.lon.max()+.10),float(sites.lat.max()+.08)]
if not clip.exists():
    h=requests.get(URL,headers={'Range':'bytes=0-16383'},timeout=30)
    h.raise_for_status()
    assert h.status_code==206 and len(h.content)==16384
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_TIMEOUT='30',CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif'):
        with rasterio.open('/vsicurl/'+URL) as src:
            q=from_bounds(*bbox,transform=src.transform)
            col,row=math.floor(q.col_off),math.floor(q.row_off)
            win=Window(col,row,math.ceil(q.col_off+q.width)-col,math.ceil(q.row_off+q.height)-row)
            data=src.read(1,window=win)
            tr=src.window_transform(win)
            profile=src.profile.copy();profile.update(width=data.shape[1],height=data.shape[0],transform=tr,nodata=255)
            with rasterio.open(clip,'w',**profile) as dst:dst.write(data,1)
            m=dict(source_url=URL,version='GSW1.5',period='1984-2024',temporal_role='longterm_background_not_2016',http_range_probe='bytes=0-16383',http_status=h.status_code,etag=h.headers.get('ETag'),last_modified=h.headers.get('Last-Modified'),remote_total_bytes=h.headers.get('Content-Range'),requested_bbox_wgs84=bbox,source_window=[win.col_off,win.row_off,win.width,win.height],native_pixel_degrees=list(src.res),nominal_resolution_m=30,source_crs=str(src.crs),source_nodata_tag=src.nodata,explicit_nodata_from_guide=255,window_sha256=hashlib.sha256(clip.read_bytes()).hexdigest())
            manifest.write_text(json.dumps(m,indent=2),encoding='utf8')
m=json.loads(manifest.read_text(encoding='utf8'))
assert hashlib.sha256(clip.read_bytes()).hexdigest()==m['window_sha256']
with rasterio.open(clip) as src:
    a=src.read(1);tr=src.transform;bounds=src.bounds
    assert set(np.unique(a)).issubset(set(range(101))|{255})
    rr,cc=np.indices(a.shape)
    lon=tr.c+(cc+.5)*tr.a;lat=tr.f+(rr+.5)*tr.e
    project=Transformer.from_crs(4326,2056,always_xy=True)
    xx,yy=project.transform(lon,lat)
    sx,sy=project.transform(sites.lon.to_numpy(),sites.lat.to_numpy())
    tree={t:cKDTree(np.column_stack([xx[(a>=t)&(a<=100)],yy[(a>=t)&(a<=100)]])) for t in [1,50]}
    # approximate local pixel area in metres; WGS84 geographic pixel widths vary with latitude.
    area=np.cos(np.deg2rad(lat))
    records=[]
    for i,(_,s) in enumerate(sites.iterrows()):
        dist=np.hypot(xx-sx[i],yy-sy[i]);search=dist<=5000
        valid=(a<=100)
        # Require full mapped search circle and no missing cells before exact nearest label.
        edge=float(min(dist[0,:].min(),dist[-1,:].min(),dist[:,0].min(),dist[:,-1].min()))
        assert edge>5000
        dx,dy=project.transform(s.lon+abs(tr.a),s.lat);width=np.hypot(dx-sx[i],dy-sy[i])
        dx,dy=project.transform(s.lon,s.lat+abs(tr.e));height=np.hypot(dx-sx[i],dy-sy[i])
        halfdiag=float(np.hypot(width,height)/2)
        r=dict(site_id=s.site_id,lat=s.lat,lon=s.lon,product='JRC_GSW_v1.5_occurrence',period='1984-2024',study_year=2016,temporal_role='background_including_post_experiment_years',annual2016_status='not_extracted_corrected_asset_requires_GEE_access',pixel_width_m=width,pixel_height_m=height,search_radius_m=5000,search_valid_fraction=float(valid[search].mean()),search_edge_distance_m=edge,coordinate_extent_scenario_m=100,pixel_halfdiagonal_m=halfdiag)
        ri,ci=src.index(s.lon,s.lat);v=int(a[ri,ci]);r['point_occurrence_percent']=v if v<=100 else np.nan
        for t in [1,50]:
            nearest=float(tree[t].query([sx[i],sy[i]])[0]);unknown=search&~valid
            missing_near=float(dist[unknown].min()) if unknown.any() else np.inf
            gate=nearest<=5000 and missing_near>nearest+halfdiag
            r[f'distance_occ_ge{t}_pixelcentre_m']=nearest if gate else np.nan
            r[f'distance_occ_ge{t}_status']='computed_detected_pixelcentre' if gate else ('right_censored_5000m' if not unknown.any() else 'missing_support')
            r[f'distance_occ_ge{t}_censor_lower_bound_m']=5000 if not gate and not unknown.any() else np.nan
            # Scenario bounds, not uncertainty confidence intervals or geolocation accuracy claims.
            r[f'distance_occ_ge{t}_scenario_min_m']=max(0,nearest-100-halfdiag) if gate else np.nan
            r[f'distance_occ_ge{t}_scenario_max_m']=nearest+100+halfdiag if gate else np.nan
        for radius in [100,500,1000]:
            mask=dist<=radius;ok=mask&valid
            coverage=float(area[ok].sum()/area[mask].sum())
            r[f'valid_area_fraction_{radius}m']=coverage
            r[f'mean_occurrence_{radius}m_percent']=float(np.average(a[ok],weights=area[ok])) if coverage>=.95 else np.nan
            for t in [1,50]:
                r[f'water_fraction_occ_ge{t}_{radius}m']=float(area[ok&(a>=t)].sum()/area[ok].sum()) if coverage>=.95 else np.nan
        records.append(r)
pd.DataFrame(records).to_csv(OUT/'giavi_water_context.csv',index=False)
(OUT/'source_manifest.json').write_text(json.dumps(m,indent=2),encoding='utf8')
print(pd.DataFrame(records)[['site_id','distance_occ_ge1_pixelcentre_m','distance_occ_ge50_pixelcentre_m','water_fraction_occ_ge1_500m','search_valid_fraction']].to_string(index=False))
