"""Synthetic GloBI complete-release schema and SHA safeguards; not ecological data."""
import sys,hashlib,tempfile,unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'/'analysis'))
import alan_globi_full_v011_extract as s

class FullStableGloBITests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        d=pd.DataFrame([
            ('Sphingidae','Lepidoptera','Macroglossum stellatarum','Animalia','Plantae','Lonicera japonica','visitsFlowersOf','field'),
            ('Nymphalidae','Lepidoptera','Vanessa cardui','Animalia','Plantae','Lonicera japonica','visitsFlowersOf','field'),
            ('Noctuidae','Lepidoptera','Autographa gamma','Animalia','Plantae','Silene latifolia','pollinates','iNaturalist'),
            ('Pieridae','Lepidoptera','Pieris rapae','Animalia','Plantae','Silene latifolia','flowersVisitedBy','field'),
        ],columns=['sourceTaxonFamilyName','sourceTaxonOrderName','sourceTaxonSpeciesName',
                   'sourceTaxonKingdomName','targetTaxonKingdomName',
                   'targetTaxonSpeciesName','interactionTypeName','studyTitle'])
        d['sourceTaxonName']=d.sourceTaxonSpeciesName
        d['targetTaxonName']=d.targetTaxonSpeciesName
        d['targetTaxonFamilyName']='Caprifoliaceae'
        d['targetTaxonOrderName']='Dipsacales'
        d['targetTaxonPathNames']='Plantae|Angiosperms'
        d['sourceTaxonPathNames']='Animalia|Lepidoptera'
        d.loc[3,'sourceTaxonFamilyName']='Caryophyllaceae'
        d.loc[3,'sourceTaxonOrderName']='Caryophyllales'
        d.loc[3,'sourceTaxonKingdomName']='Plantae'
        d.loc[3,'sourceTaxonSpeciesName']='Silene latifolia'
        d.loc[3,'sourceTaxonName']='Silene latifolia'
        d.loc[3,'sourceTaxonPathNames']='Plantae|Angiosperms'
        d.loc[3,'targetTaxonFamilyName']='Pieridae'
        d.loc[3,'targetTaxonOrderName']='Lepidoptera'
        d.loc[3,'targetTaxonKingdomName']='Animalia'
        d.loc[3,'targetTaxonSpeciesName']='Pieris rapae'
        d.loc[3,'targetTaxonName']='Pieris rapae'
        d.loc[3,'targetTaxonPathNames']='Animalia|Arthropoda|Lepidoptera'
        self.src=self.root/'interactions.csv.gz'
        d.to_csv(self.src,index=False,compression='gzip')
        self.digest=hashlib.sha256(self.src.read_bytes()).hexdigest()

    def test_reverse_and_forward_relations_both_retained(self):
        out=self.root/'out'
        report=s.extract(self.src,out,expected_sha=self.digest,allow_fixture=True,chunk_size=100)
        self.assertEqual(report['moth_qualified_rows'],2)
        self.assertEqual(report['butterfly_qualified_rows'],2)
        self.assertEqual(report['status'],'SYNTHETIC_FIXTURE_NOT_REAL')
        flags=pd.read_csv(out/'GloBI_v011_positive_species_flags.csv.gz').set_index('plant_name')
        self.assertEqual(flags.loc['Silene latifolia','all_butterfly_any_family'],1)
        self.assertEqual(flags.loc['Silene latifolia','non_inat_moth_any_family'],0)

    def test_sha_mismatch_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'SHA256 mismatch'):
            s.extract(self.src,self.root/'out',expected_sha='0'*64,allow_fixture=True)

    def test_unpinned_fixture_cannot_claim_release(self):
        with self.assertRaisesRegex(ValueError,'Cannot silently'):
            s.extract(self.src,self.root/'out',expected_sha=self.digest)

    def test_missing_family_is_not_inferred_as_moth(self):
        d=pd.read_csv(self.src)
        d.loc[0,'sourceTaxonFamilyName']=''
        d.to_csv(self.src,index=False,compression='gzip')
        sha=hashlib.sha256(self.src.read_bytes()).hexdigest()
        report=s.extract(self.src,self.root/'out',expected_sha=sha,allow_fixture=True,chunk_size=100)
        self.assertEqual(report['moth_qualified_rows'],1)

if __name__=='__main__':unittest.main()
