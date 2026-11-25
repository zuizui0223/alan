#!/usr/bin/env python3
"""Reproduce Robert order and bat worksheet extraction, independently of JS XML reader."""
import argparse, hashlib, json, zipfile
from pathlib import Path
import openpyxl
from openpyxl.utils import get_column_letter
p=argparse.ArgumentParser()
p.add_argument("xlsx",type=Path)
p.add_argument("extracted_json",type=Path)
args=p.parse_args()
b=args.xlsx.read_bytes()
assert len(b)==62795
assert hashlib.md5(b).hexdigest()=="8fd2cee16876f8fba6260aa91ea636dc"
assert hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()=="8bf0913affe686405e6edcd9696c3a6e81c8ab3f"
with zipfile.ZipFile(args.xlsx) as z:
    assert len(z.namelist())==15 and z.testzip() is None
w=openpyxl.load_workbook(args.xlsx, read_only=True,data_only=True)
stored=json.loads(args.extracted_json.read_text())
for sheet,field in [("Insects_Order","insects"),("Bat activity_Species","bat_activity")]:
    rows=[{get_column_letter(i):v for i,v in enumerate(r,1) if v is not None}
          for r in w[sheet].values]
    rows=[r for r in rows if r]
    expected=stored[field]
    assert len(rows)==len(expected),(sheet,len(rows),len(expected))
    for a,c in zip(rows,expected):
        assert a.keys()==c.keys(),(sheet,a.keys(),c.keys())
        for k,v in a.items():
            assert v==c[k] or (isinstance(v,(float,int)) and isinstance(c[k],(float,int)) and abs(v-c[k])<1e-12),(sheet,k)
print("verified both original worksheet extracts")
