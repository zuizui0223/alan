"""Compare visual transcription of original Supplement Table S1 with biomass rows.

Table cells are counts, not digitized plot estimates. No absences imputed.
"""
from pathlib import Path
import hashlib
import json
import sys
import pandas as pd

OUT = Path(__file__).resolve().parent
# Columns C:T1,T2,T3; A:T1,T2,T3; W:T1,T2,T3; printed total.
S1 = {
    'Agrostis tenuis': [6,6,6,6,6,6,6,6,6,54],
    'Anthoxanthum odoratum': [6,6,6,6,6,6,6,6,6,54],
    'Cynosurus cristatus': [0,0,1,0,1,2,0,0,2,6],
    'Dactylis glomerata': [6,6,5,6,6,6,6,6,6,53],
    'Festuca ovina': [6,6,6,4,4,6,6,6,6,50],
    'Holcus lanatus': [6,6,5,5,5,6,6,6,6,51],
    'Achillea millefolium': [6,6,6,5,6,6,6,6,6,53],
    'Centaurea nigra': [5,6,6,4,6,6,6,6,5,50],
    'Hypochaeris radicata': [6,6,6,6,6,6,6,6,6,54],
    'Leontodon saxatilis': [1,2,3,4,6,3,4,3,4,30],
    'Leucanthemum vulgare': [3,2,0,0,0,2,4,0,2,13],
    'Lotus corniculatus': [5,5,6,5,4,6,5,6,5,47],
    'Lotus pedunculatus': [6,6,6,6,6,6,6,6,6,54],
    'Plantago lanceolata': [6,5,6,5,6,6,6,6,6,52],
    'Prunella vulgaris': [5,6,6,1,6,3,3,4,5,39],
    'Ranunculus acris': [3,2,4,0,1,1,0,1,2,14],
    'Trifolium dubium': [0,0,2,1,0,0,1,0,0,4],
    'Trifolium pratense': [6,6,5,5,6,6,5,6,6,51],
}


def main(cache):
    source = Path(cache)/'anic2022_actual_supp.pdf'
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    assert digest == '9e8f61fa1f3ed396575dfaa8f19ad63c495138171fce2cdc29e678aae8a14a27'
    biomass = pd.read_csv(OUT/'species_biomass_source_rows.csv')
    rows = []
    for species, values in S1.items():
        assert sum(values[:9]) == values[9]
        for i, count in enumerate(values[:9]):
            light, trophic = f'L{i//3+1}', f'T{i%3+1}'
            excluded = species == 'Prunella vulgaris' and light == 'L2' and trophic == 'T1'
            selected = biomass[(biomass.species == species) & (biomass.light == light) & (biomass.trophic == trophic)]
            n = len(selected)
            listed = species in set(biomass.species)
            expected = count-int(excluded)
            status = ('not_in_13_species_biomass_table' if not listed else
                      'author_glm_exclusion_matches_unlisted_biomass_cell' if excluded and n == expected else
                      'agrees' if n == expected else 'unresolved_source_count_conflict')
            rows.append(dict(species=species, light=light, trophic=trophic,
                surviving_mesocosms_S1=count, printed_species_total=values[9],
                S1_excluded_from_glm=excluded, biomass_source_row_count=n,
                biomass_count_after_documented_exclusion_expected=expected,
                in_13_species_table=listed, status=status, source_location='Supplement Table S1 PDF page5',
                extraction='visual_table_transcription_not_plot_digitization'))
    result = pd.DataFrame(rows)
    result.to_csv(OUT/'species_survival_reconciliation.csv', index=False)
    conflicts = result[result.status.eq('unresolved_source_count_conflict')]
    assert len(conflicts) == 3 and conflicts.species.eq('Centaurea nigra').all()
    matched = result[result.in_13_species_table].copy()
    audit = dict(source_doi='10.1002/ece3.9157',
        source_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9352868/supplementaryFiles',
        supplement_sha256=digest,
        transcribed_cells=len(result), transcribed_species=18,
        compared_cells=len(matched), agreeing_cells=int(matched.status.eq('agrees').sum()),
        documented_exclusion_cells=1, unexplained_count_conflict_cells=len(conflicts),
        conflicting_species='Centaurea nigra: S1 total50 vs biomass54 positive rows',
        newly_estimated_survival_effects=0, imputed_biomass_values=0,
        final_adjudication=False, human_duplicate_extraction_complete=False)
    (OUT/'species_survival_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main(sys.argv[1])
