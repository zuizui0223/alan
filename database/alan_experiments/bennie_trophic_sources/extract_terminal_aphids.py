"""Audit later annual terminal aphid counts; retain undefined zero-arm ratios."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent


def main(cache):
    cache = Path(cache)
    allocation = pd.read_csv(OUT.parent/'plant_sources/anic2022_community_units.csv').set_index('mesocosm')
    final_rows, arm_rows, effect_rows, contrast_rows = [], [], [], []
    hashes = {}
    light_codes = {'Control':'L1', 'Amber':'L2', 'White':'L3'}
    trophic_codes = {'Non Grazer':'T1', 'Grazer':'T2', 'Predator':'T3'}
    for year in [2014, 2015]:
        file = cache/f'PeaAphid{year}.csv'
        hashes[file.name] = hashlib.sha256(file.read_bytes()).hexdigest()
        data = pd.read_csv(file)
        assert len(data)==54 and data.Plot.is_unique and set(data.Plot)==set(allocation.index)
        dates = pd.to_datetime(data.columns[3:], dayfirst=True)
        terminal = data.columns[3:][dates.argmax()]
        assert terminal == {2014:'29/08/2014',2015:'27/08/2015'}[year]
        for r in data.to_dict('records'):
            unit = allocation.loc[r['Plot']]
            assert unit.light==light_codes[r['L treatment']]
            assert unit.trophic==trophic_codes[r['T treatment']]
            final_rows.append(dict(year=year,date=dates.max().date().isoformat(),plot=r['Plot'],
                light=r['L treatment'],trophic=r['T treatment'],count=r[terminal],
                record_status='missing_not_zero' if pd.isna(r[terminal]) else 'observed_count'))
        arms = {}
        for (trophic,light),g in data.groupby(['T treatment','L treatment']):
            values=g[terminal].dropna()
            assert len(g)==6
            if trophic=='Non Grazer': assert len(values)==0
            else: assert len(values)==6 and values.ge(0).all()
            mean=float(values.mean()) if len(values) else np.nan
            arms[trophic,light]=mean
            arm_rows.append(dict(year=year,date=dates.max().date().isoformat(),trophic=trophic,light=light,
                allocated_mesocosms=6,observed_mesocosms=len(values),missing_mesocosms=6-len(values),
                mean=mean,total=float(values.sum()) if len(values) else np.nan,
                zero_count_mesocosms=int(values.eq(0).sum())))
        for trophic in ['Grazer','Predator']:
            for light in ['Amber','White']:
                lit,dark=arms[trophic,light],arms[trophic,'Control']
                assert dark==0
                effect_rows.append(dict(year=year,trophic=trophic,light=light,
                    mean_lit=lit,mean_dark=dark,lnRR=np.nan,sampling_variance=np.nan,
                    status='undefined_both_means_zero' if lit==0 else 'undefined_zero_dark_mean',
                    pseudocount_added=False))
        for light in ['Amber','White']:
            contrast_rows.append(dict(year=year,light=light,predator_log_response_ratio_difference=np.nan,
                sampling_variance=np.nan,status='not_estimable_from_zero_control_mean_ratios'))
    pd.DataFrame(final_rows).to_csv(OUT/'aphid_later_terminal_units.csv',index=False)
    pd.DataFrame(arm_rows).to_csv(OUT/'aphid_later_terminal_arms.csv',index=False)
    pd.DataFrame(effect_rows).to_csv(OUT/'aphid_later_terminal_effect_status.csv',index=False)
    pd.DataFrame(contrast_rows).to_csv(OUT/'aphid_later_terminal_predator_status.csv',index=False)
    audit=dict(input_sha256=hashes,years=[2014,2015],physical_mesocosms=54,
        annual_rows=108,observed_trophic_mesocosms_per_year=36,unrecorded_plants_only_per_year=18,
        finite_lnRRs=0,finite_predator_modification_contrasts=0,
        main_2013_estimate_changed=False,new_independent_programmes=0,
        scope='Final recorded annual aphid search counts, not extinction or mortality',
        temporal_selection='Existing pooled estimate uses2013 first-report terminal; later programme terminals are not estimable by lnRR',
        introduction_regime='Aphids introduced May-June2013; predators supplemented2013-2015; 2015 herbivore reintroduction in Methods refers to slugs, not aphids',
        no_pseudocounts=True,no_zero_variance_inverse_weights=True)
    (OUT/'aphid_later_terminal_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print(pd.DataFrame(arm_rows).to_string(index=False))


if __name__=='__main__':
    main(sys.argv[1])
