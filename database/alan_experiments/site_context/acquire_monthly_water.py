"""Same-month landscape context, not experimental plot wetness or water level."""
from pathlib import Path
import hashlib
import json
import math
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds, Window
from pyproj import Geod

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
CACHE = ROOT.parent / 'gbif_annual/alan_experiment_sources/monthly_water'
SCHEDULE = {'Waterfowl_Pond_Seven_Ponds': ['2022-06'],
            'Dongting_Lake_Wetland_Station': ['2023-05', '2023-06'],
            'Gansu_Gannan_station': ['2024-07', '2024-08']}
BASE = 'https://s3.waw4-1.cloudferro.com/swift/v1/global-surface-water/download2024/monthlyhistory/VER1-5'


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    sites = pd.read_csv(OUT / 'primary_coordinates.csv').set_index('site_id')
    records, manifests = [], []
    geod = Geod(ellps='WGS84')
    for site, months in SCHEDULE.items():
        s = sites.loc[site]
        lat, lon = float(s.lat), float(s.lon)
        west, north = math.floor(lon/10)*10, math.ceil(lat/10)*10
        tile = f'{abs(west)}{"E" if west >= 0 else "W"}_{abs(north)}{"N" if north >= 0 else "S"}'
        # Margin exceeds 2km at all three sites; verified geodesically below.
        bbox = [lon-.04, lat-.04, lon+.04, lat+.04]
        for month in months:
            year, mm = month.split('-')
            url = f'{BASE}/{year}/{year}_{mm}/monthlyhistory_{tile}_v1_5_{year}_{mm}.tif'
            clip = CACHE / f'{site}_{month}.tif'
            meta = clip.with_suffix('.json')
            if not clip.exists():
                with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR', GDAL_HTTP_TIMEOUT='30', CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif'):
                    with rasterio.open('/vsicurl/'+url) as src:
                        q = from_bounds(*bbox, transform=src.transform)
                        col, row = math.floor(q.col_off), math.floor(q.row_off)
                        win = Window(col, row, math.ceil(q.col_off+q.width)-col, math.ceil(q.row_off+q.height)-row)
                        assert col >= 0 and row >= 0 and col+win.width <= src.width and row+win.height <= src.height
                        data = src.read(1, window=win)
                        profile = src.profile.copy()
                        profile.update(width=data.shape[1], height=data.shape[0], transform=src.window_transform(win))
                        with rasterio.open(clip, 'w', **profile) as dst:
                            dst.write(data, 1)
                        m = dict(site_id=site, month=month, url=url, retrieved_utc=datetime.now(timezone.utc).isoformat(),
                                 bbox=bbox, window=[win.col_off, win.row_off, win.width, win.height],
                                 source_nodata=src.nodata, source_crs=str(src.crs), native_resolution=list(src.res),
                                 sha256=hashlib.sha256(clip.read_bytes()).hexdigest())
                        meta.write_text(json.dumps(m, indent=2)+'\n', encoding='utf-8')
            m = json.loads(meta.read_text(encoding='utf-8'))
            assert m['url'] == url and hashlib.sha256(clip.read_bytes()).hexdigest() == m['sha256']
            assert m['bbox'] == bbox, 'Cached extent belongs to different input coordinates'
            manifests.append(m)
            with rasterio.open(clip) as src:
                a = src.read(1)
                assert src.crs.to_epsg() == 4326
                assert set(np.unique(a)).issubset({0, 1, 2, 255})
                rr, cc = np.indices(a.shape)
                x = src.transform.c+(cc+.5)*src.transform.a
                y = src.transform.f+(rr+.5)*src.transform.e
                _, _, dist = geod.inv(np.full(x.shape, lon), np.full(y.shape, lat), x, y)
                assert min(dist[0].min(), dist[-1].min(), dist[:,0].min(), dist[:,-1].min()) > 2000
                area = np.cos(np.deg2rad(y))  # relative geographic-cell area
                for radius in [1000, 2000]:
                    mask = dist <= radius
                    valid = mask & np.isin(a, [1, 2])
                    water = mask & (a == 2)
                    coverage = float(area[valid].sum()/area[mask].sum())
                    fraction = float(area[water].sum()/area[valid].sum()) if valid.any() else None
                    records.append(dict(site_id=site, programme_id=s.programme_id, month=month,
                        radius_m=radius, total_pixels=int(mask.sum()), valid_pixels=int(valid.sum()), water_pixels=int(water.sum()),
                        valid_area_fraction=coverage, observed_water_fraction=fraction,
                        water_fraction_if_coverage_ge95=fraction if coverage >= .95 else None,
                        missingness_status='coverage_ge95' if coverage >= .95 else 'insufficient_coverage',
                        coordinate_precision=s.coordinate_precision_or_extent,
                        interpretation='representative_site_landscape_context_not_plot_hydrology'))
            print(site, month, 'acquired', flush=True)
    pd.DataFrame(records).to_csv(OUT/'monthly_water_context.csv', index=False)
    audit = dict(product='JRC_GSW_v1.5_monthlyhistory', source_guide='https://storage.googleapis.com/water-world/downloads_ancillary/DataUsersGuidev2024_v.5.pdf',
                 input_coordinates_sha256=hashlib.sha256((OUT/'primary_coordinates.csv').read_bytes()).hexdigest(),
                 site_months=len(manifests), sites=len(SCHEDULE), input_windows=manifests,
                 no_observation_codes=[0,255], valid_codes=[1,2], water_code=2,
                 limitation='Monthly detection not hydroperiod; kilometre context not plot moisture. No effect-model adjustment or additional experimental replication.')
    (OUT/'monthly_water_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
