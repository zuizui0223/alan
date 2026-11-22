"""No-network checks for the true source and independent-unit gates."""
import importlib.util
import unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'scripts/acquisition/alan_carannante2021_boris_provenance_gate.py'
s=importlib.util.spec_from_file_location('boris_source',P);module=importlib.util.module_from_spec(s);s.loader.exec_module(module)
class Tests(unittest.TestCase):
    def test_allowed_host_only(self):
        for u in module.URLS:module.check_url(u)
        for u in ('http://boris.unibe.ch/151803/','https://evil.example/','https://boris.unibe.ch.evil.example/'):
            with self.assertRaises(ValueError):module.check_url(u)
    def test_data_href_filter(self):
        h=module.Anchors();h.feed('<a href="/bitstreams/abc/download">a</a><a href="/users/1">b</a>')
        self.assertEqual(h.links,['/bitstreams/abc/download','/users/1'])
    def test_true_units(self):
        self.assertIn('4 true lit/dark',Path(P).read_text())
        self.assertIn('NOT light×distance',Path(P).read_text())
if __name__=='__main__':unittest.main()
