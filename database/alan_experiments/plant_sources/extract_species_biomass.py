"""Preserve printed species biomass rows and descriptive contrasts; no imputation."""
from pathlib import Path
import hashlib
import json
import re
import sys

import numpy as np
import pandas as pd
import pymupdf

OUT = Path(__file__).resolve().parent
EXPECTED = 'baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0'


def main(source_directory):
    source = Path(source_directory) / 'anic2022_data.pdf'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == EXPECTED
    doc = pymupdf.open(source)
    pages = [page.get_text() for page in doc]
    text = '\n'.join(pages)
    starts = np.cumsum([0] + [len(page) + 1 for page in pages])
    pattern = re.compile(
        r'([A-Z][a-z]+ [a-z]+)\s+([A-H][1-7])\s+(L[123])\s+(T[123])'
        r'\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)')
    rows = []
    for match in pattern.finditer(text):
        species, mesocosm, light, trophic, total, root, shoot = match.groups()
        rows.append(dict(species=species, mesocosm=mesocosm, light=light,
                         trophic=trophic, total_source=total, root_source=root,
                         shoot_source=shoot, total_g=float(total), root_g=float(root),
                         shoot_g=float(shoot), source_page=int(np.searchsorted(
                             starts, match.start(), side='right'))))
    data = pd.DataFrame(rows)
    assert len(data) == 665 and data.species.nunique() == 13
    assert not data.duplicated(['species', 'mesocosm']).any()
    data['component_difference_g'] = data.total_g-data.root_g-data.shoot_g
    data['component_inconsistency'] = data.component_difference_g.abs() > .015
    community = pd.read_csv(OUT/'anic2022_community_units.csv')
    joined = data.merge(community, on='mesocosm', validate='many_to_one', suffixes=('', '_community'))
    assert joined.light.eq(joined.light_community).all()
    assert joined.trophic.eq(joined.trophic_community).all()
    data.to_csv(OUT/'species_biomass_source_rows.csv', index=False)
    data[data.component_inconsistency].to_csv(OUT/'species_component_inconsistencies.csv', index=False)
    support = data.groupby('species').agg(recorded_mesocosms=('mesocosm', 'nunique'))
    support['unlisted_of_54'] = 54-support.recorded_mesocosms
    support.to_csv(OUT/'species_biomass_support.csv')
    arms = data.groupby(['species', 'light', 'trophic']).agg(
        recorded_mesocosms=('mesocosm', 'nunique'), mean_total_g=('total_g', 'mean'),
        source_conflict=('component_inconsistency', 'any')).reset_index()
    assert len(arms) == 116 and arms.recorded_mesocosms.between(3, 6).all()
    arms.to_csv(OUT/'species_biomass_arm_means.csv', index=False)
    contrasts = []
    for (species, trophic), group in arms.groupby(['species', 'trophic']):
        group = group.set_index('light')
        for light in ['L2', 'L3']:
            if light not in group.index or 'L1' not in group.index:
                contrasts.append(dict(species=species, trophic=trophic, light=light,
                    n_recorded_lit=int(group.loc[light].recorded_mesocosms) if light in group.index else 0,
                    n_recorded_dark=int(group.loc['L1'].recorded_mesocosms) if 'L1' in group.index else 0,
                    mean_lit_g=np.nan, mean_dark_g=np.nan, ratio=np.nan, log_ratio=np.nan,
                    source_conflict=False, status='unlisted_arm_not_zero_biomass_no_estimate'))
                continue
            a, b = group.loc[light], group.loc['L1']
            assert a.mean_total_g > 0 and b.mean_total_g > 0
            contrasts.append(dict(species=species, trophic=trophic, light=light,
                n_recorded_lit=int(a.recorded_mesocosms), n_recorded_dark=int(b.recorded_mesocosms),
                mean_lit_g=a.mean_total_g, mean_dark_g=b.mean_total_g,
                ratio=a.mean_total_g/b.mean_total_g,
                log_ratio=np.log(a.mean_total_g/b.mean_total_g),
                source_conflict=bool(a.source_conflict or b.source_conflict),
                status='descriptive_recorded_biomass_no_inferred_absence_no_pooling'))
    pd.DataFrame(contrasts).to_csv(OUT/'species_biomass_descriptive_contrasts.csv', index=False)
    reconciliation = community.merge(data.groupby('mesocosm').total_g.sum().rename(
        'listed_species_total_g'), on='mesocosm', validate='one_to_one')
    reconciliation['community_minus_listed_g'] = reconciliation.total_biomass_g-reconciliation.listed_species_total_g
    reconciliation.to_csv(OUT/'species_community_reconciliation.csv', index=False)
    audit = dict(source_doi='10.6084/m9.figshare.20304039', source_sha256=EXPECTED,
        source_pages=len(pages), species_rows=len(data), species=13, community_mesocosms=54,
        complete_54_record_species=int(support.recorded_mesocosms.eq(54).sum()),
        arm_cells=len(arms), planned_contrasts=len(contrasts),
        finite_descriptive_contrasts=int(pd.DataFrame(contrasts).ratio.notna().sum()),
        unlisted_arm='Prunella vulgaris, L2, T1',
        component_conflict_rows=int(data.component_inconsistency.sum()),
        largest_community_minus_listed_g=float(reconciliation.community_minus_listed_g.max()),
        imputed_absences=0, independent_programmes_added=0, pooled_effects_added=0,
        main_catalogue_changed=False, human_duplicate_extraction_complete=False,
        limits=['13 source-listed species; 5 rare species excluded from author species analysis',
                'Unlisted species-mesocosm combinations are not assigned zero',
                'Means condition on available species records; no survival or unconditional biomass effect',
                'Source total retained when total differs from root+shoot; affected contrast flagged',
                'Shared mesocosms and controls; no cross-species covariance or confidence intervals inferred',
                'Final whole-plant dry mass is not aboveground clipped biomass or pollination mediation'])
    (OUT/'species_biomass_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main(sys.argv[1])
