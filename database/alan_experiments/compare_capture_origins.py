"""Within-experiment aquatic-minus-terrestrial capture lnRR contrasts.

Not a pooled meta-analysis, population trend, or wet-versus-dry habitat test.
All available light colours and placements are retained, with shared controls.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
d = pd.read_csv(ROOT / 'effect_catalogue.csv').set_index('effect_id')
v = pd.read_csv(ROOT / 'sampling_covariance_all_extracted.csv', index_col=0)
pairs = [('Parkinson2020_littoral', 'P2020_aquatic_pan',
          'P2020_terrestrial_pan', 'freshwater_littoral', 'original_treatment')]
for habitat in ['littoral', 'riparian']:
    for colour in ['green', 'purple', 'red', 'white']:
        prefix = f'P2023_{habitat}_{colour}'
        pairs.append((prefix, prefix + '_aquatic_origin_capture',
                      prefix + '_terrestrial_origin_capture', habitat, colour))
a = np.zeros((len(pairs), len(d)))
rows = []
for i, (key, aquatic, terrestrial, habitat, colour) in enumerate(pairs):
    assert d.loc[aquatic, 'programme_id'] == d.loc[terrestrial, 'programme_id']
    a[i, d.index.get_loc(aquatic)] = 1
    a[i, d.index.get_loc(terrestrial)] = -1
    rows.append(dict(contrast_id=key, study_id=d.loc[aquatic, 'study_id'],
                     programme_id=d.loc[aquatic, 'programme_id'], habitat=habitat,
                     colour=colour, aquatic_effect_id=aquatic,
                     terrestrial_effect_id=terrestrial))
y = a @ d.lnRR.to_numpy()
c = a @ v.loc[d.index, d.index].to_numpy() @ a.T
assert np.isfinite(c).all() and np.allclose(c, c.T)
assert np.linalg.eigvalsh(c).min() > -1e-10
assert (np.diag(c) > 0).all()
out = pd.DataFrame(rows)
out['lnRR_aquatic_minus_terrestrial'] = y
out['sampling_variance'] = np.diag(c)
out['ratio_of_capture_response_ratios'] = np.exp(y)
out['normal_approx_lower'] = np.exp(y - 1.96 * np.sqrt(np.diag(c)))
out['normal_approx_upper'] = np.exp(y + 1.96 * np.sqrt(np.diag(c)))
out.to_csv(ROOT / 'capture_origin_contrasts.csv', index=False)
pd.DataFrame(c, index=out.contrast_id, columns=out.contrast_id).to_csv(
    ROOT / 'capture_origin_contrast_covariance.csv', index_label='contrast_id')
summary = dict(contrasts=len(out), programmes=out.programme_id.nunique(),
               positive_point_estimates=int((y > 0).sum()), pooled=False,
               limitations=[
                   'Nine contrasts arise from two programmes, not nine studies.',
                   'Both programmes are freshwater shoreline experiments in Michigan.',
                   'Origin is organism life-history origin, not trap habitat wetness.',
                   'Capture reflects attraction/interception, not emergence or population growth.',
                   'Delta-method normal intervals are descriptive; 2023 cells have only three traps.',
                   'Temporal windows, lamp characteristics and ecological settings differ.',
                   'All colours and placements retained; controls and endpoints are correlated.'])
(ROOT / 'capture_origin_contrast_summary.json').write_text(
    json.dumps(summary, indent=2) + '\n', encoding='utf-8')
print(out[['contrast_id', 'lnRR_aquatic_minus_terrestrial',
           'ratio_of_capture_response_ratios', 'normal_approx_lower',
           'normal_approx_upper']].to_string(index=False))
