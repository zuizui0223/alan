import unittest
import importlib.util
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'scripts/analysis/alan_pollution_meta_source_20261114.py'
s=importlib.util.spec_from_file_location('src',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class Source(unittest.TestCase):
    def test_delimiter_and_missing(self):
        fields,rows=m.parse_csv(b'study;light;yi;vi\r\na;ALAN;-.2;.1\r\nb;UV;NA;.2\r\n')
        self.assertEqual(fields,['study','light','yi','vi'])
        self.assertEqual(rows[1]['yi'],'NA')
    def test_bom_quoted_comma(self):
        fields,rows=m.parse_csv('\ufeffstudy,kind,value\n"a, b",ALAN,1\n'.encode())
        self.assertEqual(rows[0]['study'],'a, b')
    def test_source_digest(self):
        self.assertTrue(m.valid_md5(b'abc','900150983cd24fb0d6963f7d28e17f72'))
        self.assertFalse(m.valid_md5(b'abd','900150983cd24fb0d6963f7d28e17f72'))
if __name__=='__main__':unittest.main()
