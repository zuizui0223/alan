#!/usr/bin/env python3
"""Public BORIS source location audit for the paired 3/40m × LED/dark river trial.

No effects or independent-lamp counts inferred. Do not scrape restricted data,
republish sources, contact authors, or treat a PDF as a data table.
"""
import argparse
import hashlib
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse, urljoin
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

URLS = (
    'https://boris.unibe.ch/151803/',
    'https://boris-portal.unibe.ch/server/api/discover/search/objects?query=Carannante%20dataset&size=20',
    'https://boris-portal.unibe.ch/server/api/discover/search/objects?query=%222688-8319.12053%22&size=20',
)
HOSTS = {'boris.unibe.ch', 'boris-portal.unibe.ch'}
MAX_BODY = 2_000_000

class Anchors(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]
    def handle_starttag(self, tag, attrs):
        if tag.lower() == 'a':
            v=dict(attrs).get('href')
            if v:self.links.append(v)

def check_url(u):
    p=urlparse(u)
    if p.scheme!='https' or p.hostname not in HOSTS:
        raise ValueError('UNAPPROVED_SOURCE_HOST')

def probe_url(url):
    check_url(url)
    out={'url':url, 'status':'SOURCE_STOP'}
    try:
        with urlopen(Request(url,headers={'User-Agent':'ALAN-original-data-access-audit/2026-11-22'}), timeout=20) as resp:
            final=resp.geturl()
            check_url(final)
            out['final_url']=final
            out['http_status']=resp.status
            out['content_type']=resp.headers.get('content-type','')
            if int(resp.headers.get('content-length') or 0)>MAX_BODY:
                out['status']='REMOTE_PAYLOAD_TOO_LARGE';return out
            raw=resp.read(MAX_BODY+1)
            if len(raw)>MAX_BODY:
                out['status']='REMOTE_PAYLOAD_TOO_LARGE';return out
            out['bytes']=len(raw)
            out['sha256']=hashlib.sha256(raw).hexdigest()
            if 'json' in out['content_type'] or raw[:1] in (b'{',b'['):
                j=json.loads(raw)
                out['status']='PUBLIC_JSON_RECEIVED'
                out['json_top_keys']=sorted(j.keys())[:20] if isinstance(j,dict) else []
                # No blind following URLs or login/authentication requests.
                payload=json.dumps(j,ensure_ascii=False)
                out['mentions_carannante']=bool(re.search('Carannante',payload,re.I))
                out['mentions_dataset']=bool(re.search('dataset',payload,re.I))
                out['mentions_151803']=('151803' in payload)
                out['bitstreams_mentioned']=bool(re.search('bitstreams',payload,re.I))
            elif 'html' in out['content_type'] or b'<html' in raw[:300].lower():
                h=Anchors();h.feed(raw.decode('utf-8','replace'))
                linked=[]
                for a in h.links:
                    resolved=urljoin(final,a)
                    try:check_url(resolved)
                    except ValueError:continue
                    if re.search(r'data|file|bitstream|download|csv|xlsx|zip|151803',resolved,re.I):
                        linked.append(resolved)
                out['status']='PUBLIC_HTML_RECEIVED'
                out['public_source_link_candidates']=sorted(set(linked))[:60]
                out['contains_dataset_word']=bool(re.search('dataset',raw.decode('utf-8','replace'),re.I))
            else:
                out['status']='UNEXPECTED_CONTENT_TYPE_NO_DATA_ADOPTED'
            return out
    except (HTTPError, URLError, ValueError, OSError, json.JSONDecodeError) as err:
        out['stop']=('HTTP_'+str(err.code)) if isinstance(err,HTTPError) else type(err).__name__
        return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='results/alan_carannante2021_boris_source_gate.json');args=ap.parse_args()
    result={'source':'Carannante et al. 2021', 'article_doi':'10.1002/2688-8319.12053',
            'original_source_citation':'https://boris.unibe.ch/151803/',
            'true_design':'Gürbe river Study A: 4 true lit/dark plot pairs × 6 repeated nights; pairs at 3 and 40 m, two upstream/downstream blocks. Study B: only lit traps across 3–80m (NOT light×distance).',
            'endpoint':'adult aquatic insect trap capture, not emerging abundance or seed fitness',
            'original_plot_by_night_numeric_effects_extracted':False,
            'independent_programmes_added':0,
            'probes':[probe_url(u) for u in URLS]}
    p=Path(args.output);p.parent.mkdir(exist_ok=True,parents=True)
    p.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('SOURCE_SITE_STATUS',[(x['status'],x.get('stop'),x.get('final_url')) for x in result['probes']]);print('MANIFEST',p)
if __name__=='__main__':main()
