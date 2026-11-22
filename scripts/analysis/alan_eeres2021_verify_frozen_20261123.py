#!/usr/bin/env python3
"""Fail closed on original-source EERES 2021 replay against frozen 6-header-tank audit.

Only aggregate source-derived summaries are committed; original Dryad records
remain in original official archive. Whole-source SHA256 required.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

PINNED_SOURCE='496dcfa07436619084dc677317525dc674775e9965a0fdfb482cfa1cb4d8977f'
N_HEADERS=6
TOL=1e-9

def _close(a, b, key):
    if isinstance(b,bool) or isinstance(a,bool):
        if a!=b: raise AssertionError(f'{key}: boolean mismatch')
    elif isinstance(b,(int,float)):
        if not isinstance(a,(int,float)) or not math.isclose(a,b,rel_tol=TOL,abs_tol=TOL):
            raise AssertionError(f'{key}: {a} != {b}')
    elif isinstance(b,list):
        if not isinstance(a,list) or len(a)!=len(b):raise AssertionError(f'{key}: list length')
        for i,(v,w) in enumerate(zip(a,b)):_close(v,w,f'{key}[{i}]')
    elif isinstance(b,dict):
        if not isinstance(a,dict) or set(a)!=set(b):raise AssertionError(f'{key}: dictionary keys')
        for k in b:_close(a[k],b[k],f'{key}.{k}')
    elif a!=b:raise AssertionError(f'{key}: {a} != {b}')

def verify(actual, frozen):
    assert actual['status']=='SOURCE_WORKBOOK_VERIFIED'
    assert actual['original_workbook']['sha256']==PINNED_SOURCE
    assert frozen['source_sha256']==PINNED_SOURCE
    a=actual['original_workbook']['source_effects']
    assert a['status']=='SOURCE_DERIVED_HT_BLOCK_COMPARISONS'
    assert a['n_independent_header_blocks']==N_HEADERS
    assert a['joint_shared_HT']['status']=='SOURCE_DERIVED_EIGHT_JOINT_HT_CONTRAST_COVARIANCE'
    assert len(a['joint_shared_HT']['source_outcomes'])==8
    for k in ('source_grain','effect_summaries','joint_shared_HT'):
        _close(a[k],frozen[k],f'original.{k}')
    return {'status':'ALL_EERES2021_FROZEN_ORIGINAL_EFFECTS_REPRODUCED',
            'original_xlsx_sha256':PINNED_SOURCE,'n_independent_header_tanks':N_HEADERS,
            'source_effects':8,'relative_response_contrasts':8,'joint_covariance_cells':64}

def main():
    p=argparse.ArgumentParser();p.add_argument('--observed',required=True);p.add_argument('--frozen',required=True)
    args=p.parse_args()
    observed=json.loads(Path(args.observed).read_text(encoding='utf-8'))
    frozen=json.loads(Path(args.frozen).read_text(encoding='utf-8'))
    print(json.dumps(verify(observed,frozen),sort_keys=True))

if __name__=='__main__':main()
