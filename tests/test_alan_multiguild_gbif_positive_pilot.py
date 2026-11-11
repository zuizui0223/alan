"""Synthetic-only tests for multi-predator GBIF source QA, never ecology."""
import sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'/'acquisition'))
import alan_multiguild_gbif_positive_pilot as m

def fake_gbif_run(years,names,output,cap):
    name=names[0]
    output.mkdir(parents=True)
    audit=[{'year':year,'source_rows':3,'qualified_direct_records':2,
            'occupied_1degree_source_grids':1,'source_datasets':1,
            'exclusions':{'not_direct_observation':1}} for year in years]
    return {'species_audit':[{'species':name,'years':audit}]}

class Checks(unittest.TestCase):
    def test_positive_counts_do_not_imply_absence_or_local_predation(self):
        with tempfile.TemporaryDirectory() as p:
            with patch.object(m.gbif,'run',side_effect=fake_gbif_run):
                r=m.run(Path(p)/'out',years=(2019,2024),
                    selected=['Glossophaga soricina','Thereuopoda clunifera'],cap=1000)
            self.assertEqual(len(r['species_audit']),2)
            self.assertEqual(r['species_audit'][0]['years'][0]['qualified_direct_records'],2)
            self.assertFalse(r['false_absence_from_no_record'])
            self.assertFalse(r['actual_predator_on_moth_interaction_observed_at_GBIF_site'])
            self.assertTrue((Path(p)/'out'/'PUBLIC_source_support_receipt.json').exists())

    def test_failed_GBIF_not_exported_as_zero_observations(self):
        with tempfile.TemporaryDirectory() as p:
            with patch.object(m.gbif,'run',side_effect=RuntimeError('provider cap exceeded')):
                r=m.run(Path(p)/'out',selected=['Argiope bruennichi'])
            s=r['species_audit'][0]
            self.assertEqual(s['source_status'],'UNRESOLVED_source_incomplete_or_taxonomy')
            self.assertEqual(s['years'],[])
            self.assertIn('provider cap',s['error'])

    def test_no_unknown_species_or_source_overwrite(self):
        with tempfile.TemporaryDirectory() as p:
            with self.assertRaisesRegex(ValueError,'Unknown source-positive taxon'):
                m.run(Path(p)/'out',selected=['Fictionus example'])
            dest=Path(p)/'already'
            dest.mkdir()
            with self.assertRaisesRegex(ValueError,'Refuse overwrite'):
                m.run(dest,selected=['Glossophaga soricina'])

if __name__=='__main__':
    unittest.main()
