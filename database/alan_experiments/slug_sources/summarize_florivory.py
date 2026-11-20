"""Small-pair uncertainty for retained-flower damage, not reproductive success."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parent
source = ROOT/'paired_final_observed_effects.csv'
data = pd.read_csv(source)
out = data[(data.endpoint == 'flower_mean_damage_all_observed') & data.extraction.str.startswith('wild_')].copy()
assert len(out) == 2 and (out.n_pairs == 2).all()
assert np.isfinite(out[['lnRR','variance']]).all().all() and (out.variance > 0).all()
out['df'] = out.n_pairs - 1
margin = t.ppf(.975, out.df) * np.sqrt(out.variance)
out['ratio'] = np.exp(out.lnRR)
out['ci_lower'] = np.exp(out.lnRR - margin)
out['ci_upper'] = np.exp(out.lnRR + margin)
out['interval_method'] = 'paired_delta_t_n_pairs_minus1_exploratory'
out['pooled'] = False
out.to_csv(ROOT/'florivory_intervals.csv',index=False)
audit = dict(input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
             endpoint='mean damage over observed flowers, including undamaged flowers',
             rows=2, independent_programmes=1, source_holds_preserved=True,
             missing_flowers='not imputed as undamaged',
             dependence='One shared site pair across species at different dates; no zero cross-species covariance assumed',
             scope='retained sampled flowers; not whole-plant loss, attack rate, pollination or seed production',
             meta_effects_added=0)
(ROOT/'florivory_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(out[['extraction','n_pairs','ratio','ci_lower','ci_upper']].to_string(index=False))
