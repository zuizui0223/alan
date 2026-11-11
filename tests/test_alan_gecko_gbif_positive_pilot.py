"""GBIF gecko presence is positive-only, not a verified absence or predation."""
import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'/'acquisition'))
import alan_gecko_gbif_positive_pilot as g

class Response:
 def __init__(self,body):self.body=body;self.status_code=200
 def raise_for_status(self):pass
 def json(self):return self.body
class Session:
 def __init__(self):self.headers={}
 def get(self,url,params,timeout):
  if url.endswith('/species/match'):
   return Response({'canonicalName':params['name'],'rank':'SPECIES',
       'kingdom':'Animalia','status':'ACCEPTED','usageKey':1234,'confidence':99})
  year=int(params['year'])
  r={'key':year+1,'year':year,'species':'Hemidactylus frenatus',
     'basisOfRecord':'HUMAN_OBSERVATION','occurrenceStatus':'DETECTED',
     'decimalLatitude':18.2,'decimalLongitude':-77.3,
     'coordinateUncertaintyInMeters':200,'hasGeospatialIssue':False,
     'datasetKey':'mock'}
  return Response({'count':1,'results':[r] if params['offset']==0 else []})
class Tests(unittest.TestCase):
 def test_positive_source_mock(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'out'
   report=g.run([2019,2024],['Hemidactylus frenatus'],out,session=Session())
   self.assertEqual(report['documented_house_gecko_gridyear_counts'],{'2019':1,'2024':1})
   self.assertTrue(report['house_gecko_positive_grid_is_not_gecko_absence'])
   self.assertNotIn('18.2',(out/'house_gecko_positive_1degree_cells.csv').read_text())
 def test_unknown_uncertainty_is_not_precise_detection(self):
  r={'key':1,'year':2019,'species':'Hemidactylus frenatus',
     'basisOfRecord':'HUMAN_OBSERVATION','decimalLatitude':1,'decimalLongitude':2}
  self.assertEqual(g.record(r,'Hemidactylus frenatus',2019)[1],'precision_unknown')
 def test_nonhuman_or_wrong_year_rejected(self):
  r={'key':1,'year':2019,'species':'Hemidactylus frenatus',
     'basisOfRecord':'PRESERVED_SPECIMEN','decimalLatitude':1,'decimalLongitude':2,
     'coordinateUncertaintyInMeters':3}
  self.assertEqual(g.record(r,'Hemidactylus frenatus',2019)[1],'not_direct_observation')
  r['basisOfRecord']='HUMAN_OBSERVATION'
  self.assertEqual(g.record(r,'Hemidactylus frenatus',2024)[1],'not_exact_year')
 def test_error_if_gbif_query_too_large_not_truncated(self):
  class Oversized(Session):
   def get(self,url,params,timeout):
    x=super().get(url,params,timeout)
    if url.endswith('/occurrence/search'):x.body['count']=100001
    return x
  with self.assertRaisesRegex(ValueError,'never truncate'):
   g.year_source(Oversized(),1234,2019)
 def test_modern_GBIF_status_detected_and_not_detected(self):
  r={'key':1,'year':2019,'species':'Hemidactylus frenatus',
    'basisOfRecord':'HUMAN_OBSERVATION','decimalLatitude':18,
    'decimalLongitude':-77,'coordinateUncertaintyInMeters':20,
    'occurrenceStatus':'DETECTED'}
  value,why=g.record(r,'Hemidactylus frenatus',2019)
  self.assertIsNone(why)
  self.assertEqual(value[:2],(18,-77))
  r['occurrenceStatus']='NOT_DETECTED'
  self.assertEqual(g.record(r,'Hemidactylus frenatus',2019)[1],'not_present')
  r['occurrenceStatus']='PRESENT'
  self.assertIsNone(g.record(r,'Hemidactylus frenatus',2019)[1])
  r['occurrenceStatus']='ABSENT'
  self.assertEqual(g.record(r,'Hemidactylus frenatus',2019)[1],'not_present')
 def test_only_authorized_taxa(self):
  with tempfile.TemporaryDirectory() as t:
   with self.assertRaisesRegex(ValueError,'Unknown taxon'):
    g.run([2019],['Unknownus fictus'],Path(t)/'output',session=Session())
if __name__=='__main__':unittest.main()
