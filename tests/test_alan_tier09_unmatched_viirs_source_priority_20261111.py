"""Guardrails for tier09 source-only missing-VIIRS audit."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.acquisition.alan_tier09_unmatched_viirs_source_priority_20261111 import (
    require, sha, truth, unique_file,
)


class Tier09MissingVIIRSGuardrails(unittest.TestCase):
    def test_unknown_is_never_zero_or_false(self):
        with self.assertRaisesRegex(ValueError, "missing or malformed"):
            truth(pd.Series(["true", "unknown"]), "same-year VIIRS")
        self.assertEqual(truth(pd.Series(["True", "FALSE"]), "flag").tolist(),
                         [True, False])

    def test_provenance_file_must_be_unambiguous(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            with self.assertRaisesRegex(ValueError, "expected exactly one"):
                unique_file(root, "receipt.json")
            (root / "receipt.json").write_text("signed")
            self.assertEqual(unique_file(root, "receipt.json"),
                             root / "receipt.json")
            self.assertEqual(len(sha(root / "receipt.json")), 64)
            (root / "sub").mkdir()
            (root / "sub" / "receipt.json").write_text("duplicate")
            with self.assertRaisesRegex(ValueError, "expected exactly one"):
                unique_file(root, "receipt.json")

    def test_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "drift"):
            require(False, "input drift")


if __name__ == "__main__":
    unittest.main()
