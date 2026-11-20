"""Report small-sample approximate intervals for existing within-pond contrasts."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

root = Path(__file__).resolve().parent
units = pd.read_csv(root / "trap_units.csv")
existing = pd.read_csv(root / "habitat_moderation.csv")
rows = []
for record in existing.to_dict("records"):
    endpoint, colour = record["endpoint"], record["colour"]
    components, logs, sizes = [], [], []
    valid = True
    for habitat, treatment, sign in [("littoral", colour, 1), ("littoral", "dark", -1), ("riparian", colour, -1), ("riparian", "dark", 1)]:
        x = units.loc[(units.habitat == habitat) & (units.colour == treatment), endpoint].to_numpy()
        assert len(x) == 3 and np.isfinite(x).all()
        sizes.append(len(x))
        if x.mean() <= 0:
            valid = False
            continue
        logs.append(sign * np.log(x.mean()))
        components.append(x.var(ddof=1) / (len(x) * x.mean()**2))
    out = dict(record)
    out.update(ratio=np.nan, df=np.nan, ci_low=np.nan, ci_high=np.nan,
               interval_method="not_estimable_zero_mean", independent_cells=4, units_per_cell=3)
    if valid:
        estimate, variance = sum(logs), sum(components)
        assert np.isclose(estimate, record["lnRR_littoral_minus_riparian"])
        assert np.isclose(variance, record["variance"])
        df = variance**2 / sum(v*v/(n-1) for v,n in zip(components,sizes))
        half = t.ppf(.975,df)*np.sqrt(variance)
        out.update(ratio=np.exp(estimate), df=df, ci_low=np.exp(estimate-half),
                   ci_high=np.exp(estimate+half), interval_method="four_cell_delta_Welch_t_approximation_unadjusted")
    else:
        assert pd.isna(record["lnRR_littoral_minus_riparian"])
    rows.append(out)
pd.DataFrame(rows).to_csv(root / "habitat_moderation_intervals.csv", index=False)
manifest = {"inputs_sha256":{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ["trap_units.csv","habitat_moderation.csv"]},
            "contrasts":len(rows),"finite_intervals":sum(pd.notna(r["ratio"]) for r in rows),
            "independent_programmes":1,"intervals_are_approximate":True,
            "multiplicity_adjusted":False,"shared_control_covariance":"habitat_moderation_covariance.csv",
            "caveat":"n=3 per cell; delta normality and Welch df approximation uncertain; shared controls and overlapping endpoints; not wet/dry biome inference"}
(root / "habitat_interval_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
print(pd.DataFrame(rows).query("endpoint == 'all_invertebrate_capture'")[["colour","ratio","df","ci_low","ci_high"]].to_string(index=False))
