"""Accuracy against pyproj/PROJ Transverse Mercator on GRS80.

Reference easting/northing values were generated with pyproj 3.8.0
(PROJ 9.8.1, ellps=GRS80) and are embedded as constants so the tests
have no runtime dependency on pyproj or GeographicLib.
"""
import math
import unittest

import utm as UTM


# pyproj 3.8.0 / PROJ 9.8.1, +proj=utm +ellps=GRS80 (November 2026)
# (lat, lon), (easting, northing, zone), to_latlon kwargs, off_cm_deg
WIDE_ZONE_VALUES = [
    # Norway 60N 3E (zone 32, 6 deg off CM)
    (
        (60.0, 3.0),
        (165640.332103585592, 6666593.572024058551, 32),
        {'northern': True},
        6,
    ),
    # Svalbard 80N 8.999 (zone 31, ~6 deg off CM)
    (
        (80.0, 8.999),
        (616095.365436753840, 8887577.011552674696, 31),
        {'northern': True},
        6,
    ),
    # forced zone 32 at equator, 3E (6 deg off CM)
    (
        (0.0, 3.0),
        (-168881.688527512713, 0.000000000000, 32),
        {'northern': True},
        6,
    ),
    # forced zone 31 at equator, 12E (9 deg off CM)
    (
        (0.0, 12.0),
        (1505646.899515871890, 0.000000000000, 31),
        {'northern': True},
        9,
    ),
    # forced zone 32 at 60N, 0E (9 deg off CM)
    (
        (60.0, 0.0),
        (-961.404526800674, 6685590.893409115262, 32),
        {'northern': True},
        9,
    ),
    # forced zone 33 at equator, 6E (9 deg off CM)
    (
        (0.0, 6.0),
        (-505646.899515868747, 0.000000000000, 33),
        {'northern': True},
        9,
    ),
    # zone 33, 40S, 3 deg off CM (normal zone width)
    (
        (-40.0, 18.0),
        (756099.647972018458, 5567930.943215333857, 33),
        {'northern': False},
        3,
    ),
    # zone 33, 84N, on CM
    (
        (84.0, 15.0),
        (500000.000000000175, 9328093.830467360094, 33),
        {'northern': True},
        0,
    ),
    # Aachen, zone 32 (normal zone width)
    (
        (50.77534556, 6.08388667),
        (294408.662941387040, 5628897.512984829955, 32),
        {'northern': True},
        3,
    ),
]

# In a normal UTM zone (at most 3 deg from the central meridian) the
# extended series stays clearly below 1 micrometre (measured a few nm).
MAX_ERR_M_NORMAL = 1e-6
# At 6-9 deg from the central meridian the power series, even with the
# a^9/a^10 and completed a^7/a^8 terms, still truncates: measured max
# forward error about 13 um, max inverse about 0.27 mm (at 9 deg, mid
# latitudes). 0.5 mm leaves a small margin above that.
MAX_ERR_M_EDGE = 5e-4


def _max_err(off_cm):
    return MAX_ERR_M_NORMAL if off_cm <= 3 else MAX_ERR_M_EDGE


def _local_m(lat, dlat, dlon):
    return math.hypot(
        dlat * 111132.0,
        dlon * 111320.0 * math.cos(math.radians(lat)),
    )


class TestWideZoneAccuracy(unittest.TestCase):
    def test_forward(self):
        for (lat, lon), (e_ref, n_ref, zone), _, off in WIDE_ZONE_VALUES:
            easting, northing, z, _letter = UTM.from_latlon(
                lat, lon, force_zone_number=zone, strict=False
            )
            self.assertEqual(z, zone)
            err = math.hypot(easting - e_ref, northing - n_ref)
            self.assertLess(
                err, _max_err(off),
                msg=f'forward {lat}, {lon} zone {zone} off={off} err={err}',
            )

    def test_inverse(self):
        for (lat, lon), (e_ref, n_ref, zone), kw, off in WIDE_ZONE_VALUES:
            lat2, lon2 = UTM.to_latlon(
                e_ref, n_ref, zone, strict=False, **kw
            )
            err = _local_m(lat, lat2 - lat, lon2 - lon)
            self.assertLess(
                err, _max_err(off),
                msg=f'inverse {lat}, {lon} zone {zone} off={off} err={err}',
            )


if __name__ == '__main__':
    unittest.main()
