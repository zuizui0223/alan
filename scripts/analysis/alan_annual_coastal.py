"""Check adjacent coastal land cover and the one survey coordinate outside Europe."""
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import requests
from pyproj import Transformer
from alan_wet_shore import repo
out=repo/"results/all_species_forest/wet_annual"
g=pd.read_csv(out/"habitats.csv")
api="https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CLC2000_WM/MapServer/0/query"
def point(r):
    x,y=Transformer.from_crs(4326,3035,always_xy=True).transform(r.longitude,r.latitude)
    # Explicit 2 km square; all points within 1 km are included.
    box=f"{x-1000},{y-1000},{x+1000},{y+1000}"
    j=requests.get(api,params=dict(geometry=box,geometryType="esriGeometryEnvelope",inSR=3035,where="code_00 in ('421','422','423','521','522','523')",spatialRel="esriSpatialRelIntersects",returnGeometry="false",outFields="code_00",f="json"),timeout=40).json()
    return r.unit_id,";".join(sorted({f["attributes"]["code_00"] for f in j["features"]}))
with ThreadPoolExecutor(5) as pool:s=pd.DataFrame(pool.map(point,g.itertuples()),columns=["unit_id","coastal_classes_in_2km_square"])
s.to_csv(out/"coastal_screen.csv",index=False)
u="https://www.fws.gov/wetlandsmapservice/rest/services/Wetlands/MapServer/0/query"
j=requests.get(u,params=dict(geometry="-104.961468,39.732207",geometryType="esriGeometryPoint",inSR=4326,spatialRel="esriSpatialRelIntersects",outFields="*",returnGeometry="false",f="json"),timeout=25).json()
print("NWI Denver polygons",len(j["features"]),"adjacent coastal classes",s.coastal_classes_in_2km_square.value_counts().to_dict())
