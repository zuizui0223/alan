"""Annual GLAD open water in an approximately 1 km box around baseline inland-wetland points."""
from pathlib import Path
from types import ModuleType
from concurrent.futures import ThreadPoolExecutor
import subprocess
import pandas as pd
import numpy as np
from alan_wet_shore import cache,repo
out=repo/"results/all_species_forest/wet_annual"
g=pd.read_csv(out/"habitats.csv").merge(pd.DataFrame({"y":range(1999,2024)}),how="cross")
g["west"]=np.floor(g.longitude/10)*10;g["north"]=(np.floor(g.latitude/10)+1)*10
folder=cache/"freshwater_yearly";folder.mkdir(exist_ok=True)
src=ModuleType("locked_extractor")
exec(subprocess.check_output(["git","show","9b88381:scripts/analysis/alan_wet_water.py"],cwd=repo),src.__dict__)
def fetch(kv):return src.extract(*kv,folder)
with ThreadPoolExecutor(5) as pool:water=pd.concat(pool.map(fetch,g.groupby(["west","north","y"])))
assert not water.duplicated(["unit_id","y"]).any() and len(water)==len(g)
water.sort_values(["unit_id","y"]).to_csv(out/"water.csv",index=False)
