"""Retrieve public, explicitly listed EnviDat resources and record their content hashes."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlsplit
root = Path(__file__).resolve().parent
if (root/'source_manifest.json').exists():
    raise FileExistsError('Frozen source snapshot exists; copy this script to a new directory for a new retrieval.')
ids = {'566': 'alan---nature-sustainability','496': 'invertebrate-functional-groups-dataset'}
def get(item):
    package, url, filename = item
    path=root/package/filename
    path.parent.mkdir(exist_ok=True)
    try:
        with urlopen(url, timeout=40) as response:
            content = response.read()
            mime=response.headers.get('Content-Type')
            if 'text/html' in mime:
                raise ValueError('unexpected html')
        path.write_bytes(content)
        return dict(dataset=package, filename=filename, url=url, bytes=len(content), sha256=hashlib.sha256(content).hexdigest(), mime=mime)
    except Exception as error:
        return dict(dataset=package, filename=filename, url=url, error=str(error))
jobs=[]
for package, key in ids.items():
    url='https://www.envidat.ch/api/3/action/package_show?id='+key
    res=get((package,url,'package.json'))
    meta=json.loads((root/package/'package.json').read_bytes())['result']
    for row in meta['resources']:
        jobs.append((package,row['url'],urlsplit(row['url']).path.split('/')[-1]))
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
    manifest=list(pool.map(get,jobs))
(root/'source_manifest.json').write_text(json.dumps(dict(retrieval_date_jst='2026-11-22',resources=manifest), indent=2))
for res in manifest:
    print(res['dataset'],res['filename'],res.get('bytes',res.get('error')))
