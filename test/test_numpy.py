import unittest

import utm as UTM

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


class TestArrayInput(unittest.TestCase):
    def test_list_from_latlon_matches_scalar(self):
        lats = [50.77534556, 0.0, -41.28646]
        lons = [6.08388667, 3.0, 174.77623611]
        e, n, z, letters = UTM.from_latlon(lats, lons)
        self.assertEqual(len(e), 3)
        for i in range(3):
            es, ns, zs, ls = UTM.from_latlon(lats[i], lons[i])
            self.assertAlmostEqual(e[i], es, delta=1e-9)
            self.assertAlmostEqual(n[i], ns, delta=1e-9)
            self.assertEqual(z[i], zs)
            self.assertEqual(letters[i], ls)

    def test_list_to_latlon_matches_scalar(self):
        eastings = []
        northings = []
        zones = []
        letters = []
        for lat, lon in ((50.77534556, 6.08388667), (40.71435, -74.00597)):
            e, n, z, ltr = UTM.from_latlon(lat, lon)
            eastings.append(e)
            northings.append(n)
            zones.append(z)
            letters.append(ltr)
        lats, lons = UTM.to_latlon(eastings, northings, zones, letters)
        for i in range(2):
            lat, lon = UTM.to_latlon(
                eastings[i], northings[i], zones[i], letters[i]
            )
            self.assertAlmostEqual(lats[i], lat, delta=1e-12)
            self.assertAlmostEqual(lons[i], lon, delta=1e-12)

    def test_length_mismatch_raises(self):
        with self.assertRaises(ValueError):
            UTM.from_latlon([0.0, 1.0], [0.0])


@unittest.skipUnless(HAS_NUMPY, 'numpy is optional')
class TestNumPyInput(unittest.TestCase):
    def test_ndarray_roundtrip(self):
        lats = np.array([50.77534556, 0.0, -33.92486889])
        lons = np.array([6.08388667, 3.0, 18.42405500])
        e, n, z, letters = UTM.from_latlon(lats, lons)
        self.assertIsInstance(e, np.ndarray)
        self.assertIsInstance(n, np.ndarray)
        self.assertIsInstance(z, np.ndarray)
        for i in range(3):
            es, ns, zs, ls = UTM.from_latlon(float(lats[i]), float(lons[i]))
            self.assertAlmostEqual(float(e[i]), es, delta=1e-9)
            self.assertAlmostEqual(float(n[i]), ns, delta=1e-9)
            self.assertEqual(int(z[i]), zs)
            self.assertEqual(letters[i], ls)
        lats2, lons2 = UTM.to_latlon(e, n, z, letters)
        self.assertIsInstance(lats2, np.ndarray)
        self.assertTrue(np.allclose(lats2, lats, atol=1e-12))
        self.assertTrue(np.allclose(lons2, lons, atol=1e-12))

    def test_scalar_force_zone_on_array(self):
        lats = np.zeros(3)
        lons = np.array([3.0, 6.0, 9.0])
        e, n, z, _letters = UTM.from_latlon(lats, lons, force_zone_number=32)
        self.assertTrue(np.all(z == 32))
        e0, n0, z0, _ = UTM.from_latlon(0.0, 3.0, force_zone_number=32)
        self.assertAlmostEqual(float(e[0]), e0, delta=1e-9)
        self.assertEqual(int(z0), 32)


if __name__ == '__main__':
    unittest.main()
