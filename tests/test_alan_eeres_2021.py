import sys,unittest,io,zipfile,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[0]))
import alan_eeres_2021_source_gate as a

class Tests(unittest.TestCase):
 def test_expected_doi(self):self.assertEqual(a.DOI,'10.5061/dryad.rfj6q57k5')
 def test_reject_other_host(self):
  with self.assertRaisesRegex(ValueError,'UNAPPROVED_HOST'):a.grab('https://example.com/x')
 def test_relative_version_href(self):
  self.assertIn('https://datadryad.org/api/v2/versions/123',a.hrefs({'_links':{'version':{'href':'/api/v2/versions/123'}}}))
 def test_ignore_untrusted_href(self):
  self.assertNotIn('https://bad.org/api/v2/foo',a.hrefs({'href':'https://bad.org/api/v2/foo'}))
 def test_reject_wrong_excel_zip(self):
  b=io.BytesIO()
  with zipfile.ZipFile(b,'w') as z:z.writestr('notthis.csv','x')
  with self.assertRaisesRegex(ValueError,'SOURCE_WORKBOOK_NOT_UNIQUE'):a.find_archive(b.getvalue())
 def test_flume_not_treated_as_repeat(self):
  self.assertTrue('EERES' in a.run.__code__.co_consts or a.DOI.endswith('rfj6q57k5'))
if __name__=='__main__':unittest.main()

class OriginalSourceGeometry(unittest.TestCase):
 def synthetic(self):
  import openpyxl
  w=openpyxl.Workbook();dr=w.active;dr.title='drift'
  dr.append(['ID','HT','flume','side','treatment_up','treatment_down','date','time']+[f'Tax{i}' for i in range(38)])
  be=w.create_sheet('benthos')
  be.append(['ID','HT','flume','section','side','treatment']+[f'Tax{i}' for i in range(40)])
  m={'CT':1,'L':2,'F':3,'L+F':6}
  for ht in range(1,7):
   for f in range(3):
    flume=f'F{(ht-1)*3+f+1}'
    left,right=[('CT','CT'),('L+F','CT'),('L','F')][f]
    for side,tr in [('left',left),('right',right)]:
     for t in range(1,5):
      dr.append([f'{ht}{f}{side}{t}',f'ht{ht}',flume,side,tr,'CT','01/08/2021',f'T{t}',m[tr]]+[0]*37)
     for rep in range(2):
      be.append([f'{ht}{f}{side}{rep}',f'ht{ht}',flume,'up',side,'CT_up' if tr=='CT' else tr,m[tr]]+[0]*39)
    for rep in range(2):
     be.append([f'{ht}{f}down{rep}',f'ht{ht}',flume,'down','__','CT_down',1]+[0]*39)
  buf=io.BytesIO();w.save(buf);return buf.getvalue()
 def test_true_header_tank_units(self):
  res=a.source_effects(self.synthetic())
  self.assertEqual(res['status'],'SOURCE_DERIVED_HT_BLOCK_COMPARISONS')
  self.assertEqual(res['n_independent_header_blocks'],6)
  self.assertEqual(res['effect_summaries']['drift']['abund'][0]['RR'],2.0)
  self.assertEqual(res['effect_summaries']['benthos']['abund'][0]['RR'],2.0)
  self.assertEqual(res['effect_summaries']['drift']['abund'][1]['RR'],2.0)
 def test_original_rows_vs_summaries(self):
  r=a.source_effects(self.synthetic())
  self.assertEqual(r['n_independent_header_blocks'],6)
  self.assertEqual(len(r['effect_summaries']['drift']['abund']),2)
  self.assertEqual(len(r['effect_summaries']['benthos']['abund']),2)
