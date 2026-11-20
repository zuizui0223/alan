"""Audit independent lighting units and recovered artificial-prey observations."""
from pathlib import Path
import hashlib
import json
import pandas as pd

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/skyglow2023'
manifest = json.loads((OUT / 'source_manifest.json').read_text())
for item in manifest:
    raw = (CACHE / item['file']).read_bytes()
    assert len(raw) == item['bytes']
    assert 'md5:' + hashlib.md5(raw).hexdigest() == item['checksum']
d = pd.read_csv(CACHE / 'df_bites.csv')
key = ['block_ID', 'rep_ID', 'unit_ID', 'patch_ID', 'x', 'y']
assert not d.duplicated(key).any() and not d.isna().any().any()
assert d.groupby('unit_ID').treatment_lux.nunique().eq(1).all()
assert d.groupby('rep_ID').block_ID.nunique().eq(1).all()
assert d.bite_count.ge(0).all() and d.bite_count.mod(1).eq(0).all()
support = d.groupby(['unit_ID', 'treatment_lux', 'block_ID', 'rep_ID']).agg(
    recovered_dummies=('bite_count', 'size'),
    observed_bite_marks=('bite_count', 'sum')).reset_index()
assert len(support) == 48 and support.recovered_dummies.le(16).all()
support['planned_dummies'] = 16
support['unrepresented_dummies'] = 16 - support.recovered_dummies
support.to_csv(OUT / 'unit_interval_support.csv', index=False)
audit = dict(raw_rows=len(d), lighting_units=d.unit_ID.nunique(),
             dose_levels=d.treatment_lux.nunique(),
             unit_intervals=len(support), planned_dummy_exposures=48*16,
             unrepresented_dummy_exposures=int(support.unrepresented_dummies.sum()),
             missing_values_in_observed_rows=int(d.isna().sum().sum()),
             duplicate_observation_keys=0, checksum_agreement=True,
             effects_extracted=0,
             cautions=['Missing exposures are not zero bites.',
                       'Repeated intervals and dummies are not independent lighting units.',
                       'Artificial prey attacks are not adult moth mortality.',
                       'Lighting timing, indoor scope and shared programme require adjudication.'])
(OUT / 'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
print(json.dumps(audit, indent=2))
