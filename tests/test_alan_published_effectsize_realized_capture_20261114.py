import unittest
from scripts.analysis.alan_published_effectsize_realized_capture_20261114 import realized_capture_endpoint

class RealizedCaptureTests(unittest.TestCase):
    def test_handling_is_not_actual_capture(self):
        self.assertFalse(realized_capture_endpoint('Handled.prey.(No.predatitor.cue)'))
        self.assertTrue(realized_capture_endpoint('prey/min'))
        self.assertTrue(realized_capture_endpoint('fish.consumption'))
        self.assertFalse(realized_capture_endpoint('prey.attraction'))

if __name__=='__main__':unittest.main()
