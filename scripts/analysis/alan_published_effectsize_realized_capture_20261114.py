"""Separate realized capture/consumption from broader predator responsiveness.
The broader sensitivity remains in the stage-panel output for transparent audit.
"""
import argparse,json
from pathlib import Path
from scripts.analysis.alan_published_effectsize_stage_panel_20261114 import (
    read_source,source_contrast,synthesize,CAPTURE_ENDPOINTS)


def realized_capture_endpoint(name):
    return name in {'Consumption.of.insects','prey/min','fish.consumption.sandy.base',
                   'fish.consumption','small.prey.captured'}


def run(source,out):
    out=Path(out)
    if out.exists():raise ValueError('Refuse overwrite')
    rows,receipt=read_source(source)
    selected=[source_contrast(r) for r in rows if r['Delete']=='Keep' and r['Measure.description'] in CAPTURE_ENDPOINTS]
    retained=[r for r in selected if not r['flags'] and realized_capture_endpoint(r['endpoint'])]
    assert len(retained)==8 and {r['study'] for r in retained}=={'387','501'}
    result={'schema':'ALAN_realized_capture_consumption_separate_from_handling_v1',
            'source':receipt,'selected_comparison_rows':len(selected),
            'realized_after_quality_guards':len(retained),
            'n_publication_clusters':2,
            'handling_probability_excluded_from_realized_count_pool':True,
            'scenarios':[synthesize(retained,rho) for rho in (0.,.5,1.)],
            'primary_estimate_rho':.5,
            'limits':[
              'Two publication clusters, shorebird prey intake and perch amphipod consumption; not a moth-predation panel.',
              'Unknown correlations and extracted SD/N are approximations; remaining primary paired and sector covariance not reconstructed.',
              'The three-study broader index additionally includes a prey-handling probability and must not be called a common consumed-prey-count effect.',
              'No harmful global ecosystem, pollination-loss, island, or predator substitute-pollination effect inferred.']}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('REALIZED_CAPTURE_CONSUMPTION_STAGE',json.dumps(result['scenarios'][1]))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.source,a.out)
