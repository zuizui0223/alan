"""Approximate paired time difference from rounded author-reported summaries.

Conditional on harvest; not a removal probability or a seed-dispersal effect.
No raw observations or cluster-adjusted removal variances are reconstructed.
"""
from pathlib import Path
import json
import pandas as pd
from scipy.stats import t

OUT=Path(__file__).resolve().parent
reported=dict(doi='10.1111/1365-2664.12206',source_location='Results and Fig3b; dissertation Chapter2',
              mean_dark_minutes=84,mean_lit_minutes=196,sd_dark=42,sd_lit=82,
              paired_n=12,reported_signed_t_dark_minus_lit=-4.1,
              raw_data_recovered=False,conditional_on_harvest=True)
(OUT/'reported_summary.json').write_text(json.dumps(reported,indent=2)+'\n')
delta=reported['mean_lit_minutes']-reported['mean_dark_minutes']
se=abs(delta/reported['reported_signed_t_dark_minus_lit'])
half=t.ppf(.975,reported['paired_n']-1)*se
pd.DataFrame([dict(effect_id='Lewanzik2014_conditional_harvest_time',
    effect_metric='paired_mean_difference_minutes',estimate=delta,
    approximate_standard_error=se,approximate_variance=se**2,
    approximate_lower=delta-half,approximate_upper=delta+half,
    paired_n=12,df=11,conditional_on_harvest=True,
    provenance='rounded_author_summary_not_raw_reanalysis',
    main_lnRR_catalogue=False,pooled=False)]).to_csv(OUT/'conditional_timing_effect.csv',index=False)
assert delta==112 and abs(delta/se-4.1)<1e-12
print(f'Conditional lit-minus-dark delay {delta:.1f} minutes; approximate interval {delta-half:.1f} to {delta+half:.1f}.')
