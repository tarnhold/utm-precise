import os
import subprocess
import sys
import unittest


CLI = os.path.join(os.path.dirname(__file__), os.pardir, "scripts", "utm-converter")


class TestCli(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run(
            [sys.executable, CLI, *args],
            capture_output=True,
            text=True,
        )

    def test_latlon_outputs_metre_precision(self):
        result = self.run_cli("latlon", "50.77534556", "6.08388667")
        self.assertEqual(result.returncode, 0, result.stderr)
        parts = result.stdout.strip().split(",")
        self.assertEqual(len(parts), 4)
        easting, northing, zone, letter = parts
        self.assertRegex(easting, r"^\d+\.\d{6}$")
        self.assertRegex(northing, r"^\d+\.\d{6}$")
        self.assertEqual(zone, "32")
        self.assertEqual(letter, "U")
        self.assertAlmostEqual(float(easting), 294408.662941, places=6)
        self.assertAlmostEqual(float(northing), 5628897.512985, places=6)

    def test_utm_accepts_float_easting_northing(self):
        result = self.run_cli(
            "utm", "294408.6629413165", "5628897.512984848", "32", "U"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lat, lon = result.stdout.strip().split(",")
        self.assertAlmostEqual(float(lat), 50.77534556, delta=1e-9)
        self.assertAlmostEqual(float(lon), 6.08388667, delta=1e-9)

    def test_missing_subcommand_is_clean_error(self):
        result = self.run_cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        combined = result.stderr
        self.assertNotIn("NameError", combined)
        self.assertNotIn("Traceback", combined)
        self.assertIn("usage:", combined)
        self.assertTrue(
            "required" in combined.lower() or "COMMAND" in combined
        )


if __name__ == "__main__":
    unittest.main()
