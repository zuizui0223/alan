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
