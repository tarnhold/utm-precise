# cython: language_level=3, boundscheck=False, wraparound=False
# cython: cdivision=True, initializedcheck=False

from libc.math cimport sin, cos, tan, sqrt, fmod, M_PI

from cpython cimport array as carray
import array

from utm.error import OutOfRangeError

__all__ = ['to_latlon', 'from_latlon']

cdef double rad2deg = 180.0 / M_PI
cdef double deg2rad = M_PI / 180.0

cdef double K0 = 0.9996

# GRS80
cdef double R = 6378137.0
cdef double F = 1.0 / 298.257222101
# WGS84
#cdef double R = 6378137.0
#cdef double F = 1.0 / 298.257223563

# first numerical eccentricity squared (e^2)
cdef double E = (2.0 - F) * F
cdef double E2 = E ** 2
cdef double E3 = E ** 3
cdef double E4 = E ** 4
cdef double E5 = E ** 5
cdef double E6 = E ** 6
cdef double E7 = E ** 7
# second numerical eccentricity squared (e^'2)
cdef double E_P2 = E / (1.0 - E)

# third flattening n
cdef double _E = F / (2.0 - F)
cdef double _E2 = _E ** 2
cdef double _E3 = _E ** 3
cdef double _E4 = _E ** 4
cdef double _E5 = _E ** 5
cdef double _E6 = _E ** 6
cdef double _E7 = _E ** 7

# Meridian distance from latitude
#
# Deakin, R. E. (2006): Meridian Distance, School of Mathematical &
#   Geospatial Sciences, RMIT University, Melbourne, March 2006.
# Weintritt, A. (2013): So, What is Actually the Distance from the
#   Equator to the Pole? - Overview of Meridian Distance Approximations,
#   TransNav the International Journal on Marine Navigation and Safety
#   of Sea Transportation, Volume 7, Number 2, June 2013.
# Coefficients through e^14 by Taylor expansion of the meridian
# integrand and term-by-term integration (sympy).
cdef double M1 =                  (1 - E / 4 - 3 * E2 / 64 -  5 * E3 / 256 - 175 * E4 / 16384 - 441 * E5 / 65536 - 4851 * E6 / 1048576 - 14157 * E7 / 4194304)
cdef double M2 =   3. /       8 * (    E +         E2 /  4 + 15 * E3 / 128 +  35 * E4 /   512 + 735 * E5 / 16384 + 2079 * E6 / 65536 + 99099 * E7 / 4194304)
cdef double M3 =  15. /     256 * (                E2      +  3 * E3 /   4 +  35 * E4 /    64 + 105 * E5 /   256 + 10395 * E6 / 32768 + 33033 * E7 / 131072)
cdef double M4 =  35. /    3072 * (                               E3       +   5 * E4 /     4 + 315 * E5 /   256 +  1155 * E6 /  1024 + 33033 * E7 / 32768)
cdef double M5 = 315. /  131072 * (                                                E4         +   7 * E5 /     4 +   693 * E6 /   320 +  3003 * E7 /  1280)
cdef double M6 = 693. / 1310720 * (                                                               E5              +     9 * E6 /     4 +   429 * E7 /   128)
cdef double M7 = 1001. / 8388608 * (                                                                              E6                     +    11 * E7 /     4)
cdef double M8 = 6435. / 234881024 *                                                                                                        E7

# Latitude from meridian distance
#
# Deakin, R. E. (2012): Great Elliptic Arc Distance, School of Mathematical &
#   Geospatial Sciences, RMIT University, Melbourne, January 2012.
# Coefficients through n^7 by Fourier reversion of the meridian
# series (sympy), reducing to Deakin through n^5.
cdef double P2 = (   3. /     2 * _E  -  27. /    32 * _E3 + 269. /   512 * _E5 - 6607. / 24576 * _E7)
cdef double P3 = (  21. /    16 * _E2 -  55. /    32 * _E4 + 6759. /  4096 * _E6)
cdef double P4 = ( 151. /    96 * _E3 - 417. /   128 * _E5 + 87963. / 20480 * _E7)
cdef double P5 = (1097. /   512 * _E4 - 15543. / 2560 * _E6)
cdef double P6 = (8011. /  2560 * _E5 - 69119. / 6144 * _E7)
cdef double P7 = (293393. / 61440 * _E6)
cdef double P8 = (6459601. / 860160 * _E7)

ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"


cdef inline double mod_angle_c(double value) noexcept:
    """Angle in radians in range [-pi, pi)."""
    cdef double rem = fmod(value + M_PI, 2.0 * M_PI)
    if rem < 0.0:
        rem += 2.0 * M_PI
    return rem - M_PI


cpdef mod_angle(double value):
    """Returns angle in radians in range [-pi, pi)"""
    return mod_angle_c(value)


cdef void tm_forward(double lat_rad, double dlon, double *x, double *y) noexcept:
    """Transverse Mercator forward (metres from CM / equator, no false offsets)."""
    cdef double lat_tan, lat_tan2, lat_tan4, lat_tan6, lat_tan8
    cdef double n, c, c2, c3, c4, c5, c6, c7, c8
    cdef double a, a2, a3, a4, a5, a6, a7, a8, a9, a10
    cdef double m

    lat_tan = tan(lat_rad)
    lat_tan2 = lat_tan ** 2
    lat_tan4 = lat_tan ** 4
    lat_tan6 = lat_tan ** 6
    lat_tan8 = lat_tan ** 8

    # Normalkruemmungsradius
    n = R / sqrt(1 - E * sin(lat_rad) ** 2)
    c = E_P2 * cos(lat_rad) ** 2
    c2 = c ** 2
    c3 = c ** 3
    c4 = c ** 4
    c5 = c ** 5
    c6 = c ** 6
    c7 = c ** 7
    c8 = c ** 8

    a = cos(lat_rad) * dlon
    a2 = a ** 2
    a3 = a ** 3
    a4 = a ** 4
    a5 = a ** 5
    a6 = a ** 6
    a7 = a ** 7
    a8 = a ** 8
    a9 = a ** 9
    a10 = a ** 10

    # meridian distance from latitude
    m = R * (M1 * lat_rad -
             M2 * sin(2 * lat_rad) +
             M3 * sin(4 * lat_rad) -
             M4 * sin(6 * lat_rad) +
             M5 * sin(8 * lat_rad) -
             M6 * sin(10 * lat_rad) +
             M7 * sin(12 * lat_rad) -
             M8 * sin(14 * lat_rad))

    # Kelly, Kevin M. (1986): Coordinate Transformations - Universal
    #    Transverse Mercator/Geographic. Ontario Ministry of Natural
    #    Resources. March 1986.
    # Snyder, John P. (1987): Map Projections - A Working Manual, U.S.
    #   Geological Survey Professional Paper 1395, p.60ff, Washington.
    # Hofmann-Wellenhof, B.; Kienast, G.; Lichtenegger, H. (1994): GPS in der
    #    Praxis, Springer-Verlag Wien New York, p.97ff, Wien.
    # Terms through a^10, and the O(e'^2) parts of a^7 and a^8 omitted
    # in those sources, by Taylor expansion of the complex meridian
    # mapping F(psi + i*lambda) (sympy).
    x[0] = K0 * n * (
                        a +
                        a3 / 6 * (1 - lat_tan2 + c) +
                        a5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 14 * c + 13 * c2 + 4 * c3 - 58 * c * lat_tan2 - 64 * c2 * lat_tan2 - 24 * lat_tan2 * c3) +
                        a7 / 5040 * (61 - 479 * lat_tan2 + 179 * lat_tan4 - lat_tan6 + 331 * c + 715 * c2 + 769 * c3 + 412 * c4 + 88 * c5 - 3298 * c * lat_tan2 - 8655 * c2 * lat_tan2 - 10964 * c3 * lat_tan2 - 6760 * c4 * lat_tan2 - 1632 * c5 * lat_tan2 + 1771 * c * lat_tan4 + 6080 * c2 * lat_tan4 + 9480 * c3 * lat_tan4 + 6912 * c4 * lat_tan4 + 1920 * c5 * lat_tan4) +
                        a9 / 362880 * (1385 - 19028 * lat_tan2 + 18270 * lat_tan4 - 1636 * lat_tan6 + lat_tan8 + 12284 * c + 45318 * c2 + 90804 * c3 + 107241 * c4 + 74904 * c5 + 28720 * c6 + 4672 * c7 - 214140 * c * lat_tan2 - 951468 * c2 * lat_tan2 - 2220708 * c3 * lat_tan2 - 2984976 * c4 * lat_tan2 - 2332848 * c5 * lat_tan2 - 987904 * c6 * lat_tan2 - 175680 * c7 * lat_tan2 + 290868 * c * lat_tan4 + 1652910 * c2 * lat_tan4 + 4662840 * c3 * lat_tan4 + 7309344 * c4 * lat_tan4 + 6500736 * c5 * lat_tan4 + 3077568 * c6 * lat_tan4 + 603648 * c7 * lat_tan4 - 47188 * c * lat_tan6 - 384384 * c2 * lat_tan6 - 1394064 * c3 * lat_tan6 - 2644992 * c4 * lat_tan6 - 2741760 * c5 * lat_tan6 - 1474560 * c6 * lat_tan6 - 322560 * c7 * lat_tan6)
                        )

    y[0] = K0 * (m + n * lat_tan * (
                                        a2 / 2 +
                                        a4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c2) +
                                        a6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 270 * c + 445 * c2 + 324 * c3 + 88 * c4 - 330 * c * lat_tan2 - 680 * lat_tan2 * c2 - 600 * lat_tan2 * c3 - 192 * lat_tan2 * c4) +
                                        a8 / 40320 * (1385 - 3111 * lat_tan2 + 543 * lat_tan4 - lat_tan6 + 10899 * c + 34419 * c2 + 56385 * c3 + 50856 * c4 + 24048 * c5 + 4672 * c6 - 32802 * c * lat_tan2 - 129087 * c2 * lat_tan2 - 252084 * c3 * lat_tan2 - 263088 * c4 * lat_tan2 - 140928 * c5 * lat_tan2 - 30528 * c6 * lat_tan2 + 9219 * c * lat_tan4 + 49644 * c2 * lat_tan4 + 121800 * c3 * lat_tan4 + 151872 * c4 * lat_tan4 + 94080 * c5 * lat_tan4 + 23040 * c6 * lat_tan4) +
                                        a10 / 3628800 * (50521 - 206276 * lat_tan2 + 101166 * lat_tan4 - 4916 * lat_tan6 + lat_tan8 + 612540 * c + 3043190 * c2 + 8250228 * c3 + 13505721 * c4 + 13774680 * c5 + 8592896 * c6 + 3009024 * c7 + 454144 * c8 - 3277980 * c * lat_tan2 - 19954380 * c2 * lat_tan2 - 63646884 * c3 * lat_tan2 - 119311872 * c4 * lat_tan2 - 136688400 * c5 * lat_tan2 - 94408960 * c6 * lat_tan2 - 36196416 * c7 * lat_tan2 - 5928192 * c8 * lat_tan2 + 2402100 * c * lat_tan4 + 19210830 * c2 * lat_tan4 + 75185208 * c3 * lat_tan4 + 165884208 * c4 * lat_tan4 + 217543680 * c5 * lat_tan4 + 168626880 * c6 * lat_tan4 + 71493120 * c7 * lat_tan4 + 12801024 * c8 * lat_tan4 - 239220 * c * lat_tan6 - 2879440 * c2 * lat_tan6 - 14852880 * c3 * lat_tan6 - 40247424 * c4 * lat_tan6 - 62092800 * c5 * lat_tan6 - 55019520 * c6 * lat_tan6 - 26127360 * c7 * lat_tan6 - 5160960 * c8 * lat_tan6)
                                        ))


cdef void tm_reverse(double x, double y, double *lat_rad, double *dlon) noexcept:
    """Transverse Mercator inverse (metres from CM / equator, no false offsets)."""
    cdef double m, mu, p_rad
    cdef double p_sin2, p_cos, p_tan, p_tan2, p_tan4, p_tan6, p_tan8
    cdef double n, r, c, c2, c3, c4, c5, c6, c7, c8
    cdef double d, d2, d3, d4, d5, d6, d7, d8, d9, d10

    m = y / K0
    mu = m / (R * M1)

    # meridian distance to latitude
    p_rad = (mu +
             P2 * sin(2 * mu) +
             P3 * sin(4 * mu) +
             P4 * sin(6 * mu) +
             P5 * sin(8 * mu) +
             P6 * sin(10 * mu) +
             P7 * sin(12 * mu) +
             P8 * sin(14 * mu))

    p_sin2 = sin(p_rad) ** 2
    p_cos = cos(p_rad)

    p_tan = tan(p_rad)
    p_tan2 = p_tan ** 2
    p_tan4 = p_tan ** 4
    p_tan6 = p_tan ** 6
    p_tan8 = p_tan ** 8

    n = R / sqrt(1 - E * p_sin2)
    r = (1 - E) / (1 - E * p_sin2)

    c = E_P2 * p_cos ** 2
    c2 = c ** 2
    c3 = c ** 3
    c4 = c ** 4
    c5 = c ** 5
    c6 = c ** 6
    c7 = c ** 7
    c8 = c ** 8

    d = x / (n * K0)
    d2 = d ** 2
    d3 = d ** 3
    d4 = d ** 4
    d5 = d ** 5
    d6 = d ** 6
    d7 = d ** 7
    d8 = d ** 8
    d9 = d ** 9
    d10 = d ** 10

    # Kelly, Kevin M. (1986): Coordinate Transformations - Universal
    #   Transverse Mercator/Geographic. Ontario Ministry of Natural
    #   Resources. March 1986.
    # Snyder, John P. (1987): Map Projections - A Working Manual, U.S.
    #   Geological Survey Professional Paper 1395, p.60ff, Washington.
    # Hofmann-Wellenhof, B.; Kienast, G.; Lichtenegger, H. (1994): GPS in der
    #   Praxis, Springer-Verlag Wien New York, p.97ff, Wien.
    # Inverse series through d^10 by reversion of the forward mapping
    # (sympy), including the O(e'^2) parts of d^7 and d^8 omitted in
    # those sources.
    lat_rad[0] = (p_rad - (p_tan / r) *
                (d2 / 2 -
                 d4 / 24 * (5 + 3 * p_tan2 + c - 4 * c2 - 9 * p_tan2 * c) +
                 d6 / 720 * (61 + 90 * p_tan2 + 45 * p_tan4 + 46 * c - 3 * c2 + 100 * c3 + 88 * c4 - 252 * p_tan2 * c - 66 * p_tan2 * c2 + 84 * p_tan2 * c3 - 192 * p_tan2 * c4 - 90 * p_tan4 * c + 225 * p_tan4 * c2) -
                 d8 / 40320 * (1385 + 3633 * p_tan2 + 4095 * p_tan4 + 1575 * p_tan6 + 1731 * c - 573 * c2 - 2927 * c3 - 8808 * c4 - 11472 * c5 - 4672 * c6 - 9381 * p_tan2 * c - 8445 * p_tan2 * c2 + 9609 * p_tan2 * c3 + 19392 * p_tan2 * c4 + 44880 * p_tan2 * c5 + 30528 * p_tan2 * c6 - 7371 * p_tan4 * c + 20325 * p_tan4 * c2 + 10551 * p_tan4 * c3 - 7032 * p_tan4 * c4 - 8832 * p_tan4 * c5 - 23040 * p_tan4 * c6 - 2835 * p_tan6 * c + 4725 * p_tan6 * c2 - 11025 * p_tan6 * c3) +
                 d10 / 3628800 * (50521 + 204180 * p_tan2 + 383670 * p_tan4 + 321300 * p_tan6 + 99225 * p_tan8 + 88412 * c + 8950 * c2 + 47924 * c3 + 556345 * c4 + 1659896 * c5 + 2428672 * c6 + 1702400 * c7 + 454144 * c8 - 468864 * p_tan2 * c - 796176 * p_tan2 * c2 + 233688 * p_tan2 * c3 - 1648260 * p_tan2 * c4 - 8916120 * p_tan2 * c5 - 18142464 * p_tan2 * c6 - 17159616 * p_tan2 * c7 - 5928192 * p_tan2 * c8 - 595080 * p_tan4 * c + 1769004 * p_tan4 * c2 + 1726704 * p_tan4 * c3 - 1743570 * p_tan4 * c4 + 3306600 * p_tan4 * c5 + 16450560 * p_tan4 * c6 + 25222464 * p_tan4 * c7 + 12801024 * p_tan4 * c8 - 514080 * p_tan6 * c + 868320 * p_tan6 * c2 - 2228760 * p_tan6 * c3 - 1065780 * p_tan6 * c4 + 1242360 * p_tan6 * c5 - 1739520 * p_tan6 * c6 - 5276160 * p_tan6 * c7 - 5160960 * p_tan6 * c8 - 170100 * p_tan8 * c + 255150 * p_tan8 * c2 - 396900 * p_tan8 * c3 + 893025 * p_tan8 * c4)))

    dlon[0] = (d -
                 d3 / 6 * (1 + 2 * p_tan2 + c) +
                 d5 / 120 * (5 + 28 * p_tan2 + 24 * p_tan4 + 6 * c - 3 * c2 - 4 * c3 + 8 * p_tan2 * c + 4 * p_tan2 * c2 + 24 * p_tan2 * c3) -
                 d7 / 5040 * (61 + 662 * p_tan2 + 1320 * p_tan4 + 720 * p_tan6 + 107 * c + 43 * c2 + 97 * c3 + 188 * c4 + 88 * c5 + 440 * p_tan2 * c - 234 * p_tan2 * c2 - 772 * p_tan2 * c3 - 2392 * p_tan2 * c4 - 1632 * p_tan2 * c5 + 336 * p_tan4 * c - 192 * p_tan4 * c2 + 408 * p_tan4 * c3 + 1536 * p_tan4 * c4 + 1920 * p_tan4 * c5) +
                 d9 / 362880 * (1385 + 24568 * p_tan2 + 83664 * p_tan4 + 100800 * p_tan6 + 40320 * p_tan8 + 3116 * c + 1158 * c2 - 3500 * c3 - 11735 * c4 - 20280 * c5 - 16144 * c6 - 4672 * c7 + 26736 * p_tan2 * c - 4944 * p_tan2 * c2 + 27104 * p_tan2 * c3 + 185496 * p_tan2 * c4 + 459312 * p_tan2 * c5 + 483712 * p_tan2 * c6 + 175680 * p_tan2 * c7 + 47808 * p_tan4 * c - 20736 * p_tan4 * c2 + 576 * p_tan4 * c3 - 195984 * p_tan4 * c4 - 817344 * p_tan4 * c5 - 1239552 * p_tan4 * c6 - 603648 * p_tan4 * c7 + 24192 * p_tan6 * c - 13824 * p_tan6 * c2 + 9216 * p_tan6 * c3 + 9792 * p_tan6 * c4 + 178560 * p_tan6 * c5 + 437760 * p_tan6 * c6 + 322560 * p_tan6 * c7)) / p_cos


def latitude_to_zone_letter(double latitude):
    if -80.0 <= latitude <= 84.0:
        return ZONE_LETTERS[int(latitude + 80.0) >> 3]
    else:
        return None


def latlon_to_zone_number(double latitude, double longitude):
    if 56.0 <= latitude < 64.0 and 3.0 <= longitude < 12.0:
        return 32

    if 72.0 <= latitude <= 84.0 and longitude >= 0.0:
        if longitude < 9.0:
            return 31
        elif longitude < 21.0:
            return 33
        elif longitude < 33.0:
            return 35
        elif longitude < 42.0:
            return 37

    return int((longitude + 180.0) / 6.0) + 1


def zone_number_to_central_longitude(int zone_number):
    return (zone_number - 1) * 6 - 180 + 3


cdef inline bint _is_scalar(object value):
    if isinstance(value, (bytes, str)):
        return False
    shape = getattr(value, "shape", None)
    if shape is not None:
        return shape == ()
    try:
        len(value)
    except TypeError:
        return True
    return False


cdef object _maybe_numpy_array(object sample, object seq, object dtype=None):
    if getattr(sample, "ndim", None) is None:
        return list(seq)
    try:
        import numpy as np
    except ImportError:
        return list(seq)
    if dtype is None:
        return np.asarray(seq)
    return np.asarray(seq, dtype=dtype)


def from_latlon(latitude, longitude, force_zone_number=None, strict=True):
    """This function convert Latitude and Longitude to UTM coordinate

        Parameters
        ----------
        latitude: float or array-like
            Latitude between 80 deg S and 84 deg N, e.g. [-80.0 to 84.0]

        longitude: float or array-like
            Longitude between 180 deg W and 180 deg E, e.g. [-180.0 to 180.0).
            Array-like inputs (lists or NumPy arrays) are converted element-wise.
            NumPy is optional; arrays work without it and return lists, or
            ndarrays if a NumPy array was passed in.

        force_zone number: int or array-like
            Zone Number is represented with global map numbers of an UTM Zone
            Numbers Map. You may force conversion including one UTM Zone Number.
            More information see utmzones [1]_

        strict: bool
            Check for valid latitude and longitude value. Default is True.

       .. _[1]: http://www.jaworski.ca/utmzones.htm
    """
    # Fast path for the usual scalar call (avoid array-detection overhead).
    if isinstance(latitude, (int, float)) and isinstance(longitude, (int, float)):
        return _from_latlon_scalar(latitude, longitude, force_zone_number, strict)
    if _is_scalar(latitude) and _is_scalar(longitude):
        return _from_latlon_scalar(float(latitude), float(longitude),
                                   force_zone_number, strict)
    return _from_latlon_array(latitude, longitude, force_zone_number, strict)


cdef tuple _from_latlon_scalar(double latitude, double longitude,
                               object force_zone_number, bint strict):
    cdef int zone_number
    cdef double x, y, dlon
    cdef object zone_letter

    if strict:
        if not -80.0 <= latitude <= 84.0:
            raise OutOfRangeError('latitude out of range (must be between 80 deg S and 84 deg N)')
    # in any case we have to be strict with longitude, because we calculate
    # zone_number from it
    if not -180.0 <= longitude < 180.0:
        raise OutOfRangeError('longitude out of range (must be between 180 deg W and 180 deg E)')

    if force_zone_number is None:
        zone_number = latlon_to_zone_number(latitude, longitude)
    else:
        zone_number = force_zone_number

    if not 1 <= zone_number <= 60:
        raise OutOfRangeError('zone number out of range (must be between 1 and 60)')

    zone_letter = latitude_to_zone_letter(latitude)

    dlon = mod_angle_c(deg2rad * (longitude - zone_number_to_central_longitude(zone_number)))
    tm_forward(deg2rad * latitude, dlon, &x, &y)
    x += 500000.0
    if latitude < 0.0:
        y += 10000000.0
    return x, y, zone_number, zone_letter


def _from_latlon_array(latitude, longitude, force_zone_number, strict):
    cdef Py_ssize_t i, n
    cdef double lat, lon, x, y, dlon
    cdef int zone_number
    cdef bint has_force_arr, has_force_s
    cdef int force_s
    cdef carray.array lat_arr, lon_arr, e_arr, n_arr, z_arr
    cdef double[::1] lat_v, lon_v, e_v, n_v
    cdef long[::1] z_v

    if _is_scalar(latitude):
        latitude = [latitude]
    if _is_scalar(longitude):
        longitude = [longitude]

    try:
        n = len(latitude)
        if len(longitude) != n:
            raise ValueError('latitude and longitude must have the same length')
    except TypeError:
        raise TypeError('latitude and longitude must be scalars or array-like')

    lat_arr = array.array('d', latitude)
    lon_arr = array.array('d', longitude)
    e_arr = array.array('d', [0.0]) * n
    n_arr = array.array('d', [0.0]) * n
    z_arr = array.array('l', [0]) * n
    letters = [None] * n

    lat_v = lat_arr
    lon_v = lon_arr
    e_v = e_arr
    n_v = n_arr
    z_v = z_arr

    has_force_arr = False
    has_force_s = False
    force_s = 0
    force_zones = None
    if force_zone_number is not None:
        if _is_scalar(force_zone_number):
            has_force_s = True
            force_s = force_zone_number
            if not 1 <= force_s <= 60:
                raise OutOfRangeError(
                    'zone number out of range (must be between 1 and 60)'
                )
        else:
            if len(force_zone_number) != n:
                raise ValueError('force_zone_number must have the same length')
            has_force_arr = True
            force_zones = array.array('l', force_zone_number)

    for i in range(n):
        lat = lat_v[i]
        lon = lon_v[i]
        if strict:
            if not -80.0 <= lat <= 84.0:
                raise OutOfRangeError(
                    'latitude out of range (must be between 80 deg S and 84 deg N)'
                )
        if not -180.0 <= lon < 180.0:
            raise OutOfRangeError(
                'longitude out of range (must be between 180 deg W and 180 deg E)'
            )
        if has_force_s:
            zone_number = force_s
        elif has_force_arr:
            zone_number = force_zones[i]
            if not 1 <= zone_number <= 60:
                raise OutOfRangeError(
                    'zone number out of range (must be between 1 and 60)'
                )
        else:
            zone_number = latlon_to_zone_number(lat, lon)
            if not 1 <= zone_number <= 60:
                raise OutOfRangeError(
                    'zone number out of range (must be between 1 and 60)'
                )
        dlon = mod_angle_c(deg2rad * (lon - zone_number_to_central_longitude(zone_number)))
        tm_forward(deg2rad * lat, dlon, &x, &y)
        x += 500000.0
        if lat < 0.0:
            y += 10000000.0
        e_v[i] = x
        n_v[i] = y
        z_v[i] = zone_number
        letters[i] = latitude_to_zone_letter(lat)

    sample = latitude if getattr(latitude, "ndim", None) is not None else longitude
    return (
        _maybe_numpy_array(sample, e_arr),
        _maybe_numpy_array(sample, n_arr),
        _maybe_numpy_array(sample, z_arr),
        _maybe_numpy_array(sample, letters, dtype=object),
    )


def to_latlon(easting, northing, zone_number, zone_letter=None, northern=None,
              strict=True):
    """This function convert an UTM coordinate into Latitude and Longitude

        Parameters
        ----------
        easting: int or array-like
            Easting value of UTM coordinate

        northing: int or array-like
            Northing value of UTM coordinate
            Array-like inputs are converted element-wise. NumPy is optional.

        zone number: int or array-like
            Zone Number is represented with global map numbers of an UTM Zone
            Numbers Map. More information see utmzones [1]_

        zone_letter: str
            Zone Letter can be represented as string values. Where UTM Zone
            Designators can be accessed in [1]_

        northern: bool
            You can set True or False to set this parameter. Default is None

        strict: bool
            Check for valid easting and northing value. Default is True.

       .. _[1]: http://www.jaworski.ca/utmzones.htm

    """
    if zone_letter is None:
        letter_missing = True
    elif isinstance(zone_letter, str):
        letter_missing = (zone_letter == '')
    else:
        letter_missing = False
    if letter_missing and northern is None:
        raise ValueError('either zone_letter or northern needs to be set')
    if (not letter_missing) and northern is not None and (
            isinstance(zone_letter, str) or _is_scalar(zone_letter)):
        raise ValueError('set either zone_letter or northern, but not both')

    if (isinstance(easting, (int, float))
            and isinstance(northing, (int, float))
            and isinstance(zone_number, (int, float))):
        return _to_latlon_scalar(
            easting, northing, zone_number, zone_letter, northern, strict
        )
    if (_is_scalar(easting) and _is_scalar(northing)
            and _is_scalar(zone_number)):
        return _to_latlon_scalar(
            float(easting), float(northing), int(zone_number),
            zone_letter, northern, strict
        )
    return _to_latlon_array(
        easting, northing, zone_number, zone_letter, northern, strict
    )


cdef tuple _to_latlon_scalar(double easting, double northing, int zone_number,
                             object zone_letter, object northern, bint strict):
    cdef double x, y, lat, dlon

    north = _hemisphere(zone_letter, northern, zone_number)
    if strict:
        if not 100000.0 <= easting < 1000000.0:
            raise OutOfRangeError('easting out of range (must be between 100.000 m and 999.999 m)')
        if not 0.0 <= northing <= 10000000.0:
            raise OutOfRangeError('northing out of range (must be between 0 m and 10.000.000 m)')
    if not 1 <= zone_number <= 60:
        raise OutOfRangeError('zone number out of range (must be between 1 and 60)')

    x = easting - 500000.0
    y = northing
    if not north:
        y -= 10000000.0
    tm_reverse(x, y, &lat, &dlon)
    dlon = mod_angle_c(dlon + deg2rad * zone_number_to_central_longitude(zone_number))
    return (rad2deg * lat, rad2deg * dlon)


cdef bint _hemisphere(object zone_letter, object northern, int zone_number):
    if zone_letter is not None and zone_letter != '':
        zone_letter = zone_letter.upper()
        if len(zone_letter) > 1:
            raise OutOfRangeError('zone letter out of range (must be one character)')
        if not 'C' <= zone_letter <= 'X' or zone_letter in ['I', 'O']:
            raise OutOfRangeError('zone letter out of range (must be between C and X)')
        if zone_letter == 'X' and zone_number in [32, 34, 36]:
            raise OutOfRangeError('zone number not valid (there is no zone ' + str(zone_number) + ' in X)')
        return zone_letter >= 'N'
    return bool(northern)


def _to_latlon_array(easting, northing, zone_number, zone_letter, northern, strict):
    cdef Py_ssize_t i, n
    cdef double e, nn, x, y, lat, dlon
    cdef int zone
    cdef bint north, has_zone_arr
    cdef carray.array e_arr, n_arr, lat_arr, lon_arr
    cdef double[::1] e_v, n_v, lat_v, lon_v

    if _is_scalar(easting):
        easting = [easting]
    if _is_scalar(northing):
        northing = [northing]
    n = len(easting)
    if len(northing) != n:
        raise ValueError('easting and northing must have the same length')

    e_arr = array.array('d', easting)
    n_arr = array.array('d', northing)
    lat_arr = array.array('d', [0.0]) * n
    lon_arr = array.array('d', [0.0]) * n
    e_v = e_arr
    n_v = n_arr
    lat_v = lat_arr
    lon_v = lon_arr

    has_zone_arr = not _is_scalar(zone_number)
    if has_zone_arr and len(zone_number) != n:
        raise ValueError('zone_number must have the same length')

    # Hemisphere from a scalar zone letter / northern flag (usual case).
    # Per-point letters are accepted as a sequence of the same length.
    letter_arr = None
    if (zone_letter is not None
            and not isinstance(zone_letter, str)
            and not _is_scalar(zone_letter)):
        if len(zone_letter) != n:
            raise ValueError('zone_letter must have the same length')
        letter_arr = list(zone_letter)

    for i in range(n):
        e = e_v[i]
        nn = n_v[i]
        zone = zone_number[i] if has_zone_arr else zone_number
        if letter_arr is not None:
            north = _hemisphere(letter_arr[i], None, zone)
        else:
            north = _hemisphere(zone_letter, northern, zone)
        if strict:
            if not 100000.0 <= e < 1000000.0:
                raise OutOfRangeError(
                    'easting out of range (must be between 100.000 m and 999.999 m)'
                )
            if not 0.0 <= nn <= 10000000.0:
                raise OutOfRangeError(
                    'northing out of range (must be between 0 m and 10.000.000 m)'
                )
        if not 1 <= zone <= 60:
            raise OutOfRangeError(
                'zone number out of range (must be between 1 and 60)'
            )
        x = e - 500000.0
        y = nn
        if not north:
            y -= 10000000.0
        tm_reverse(x, y, &lat, &dlon)
        lat_v[i] = rad2deg * lat
        lon_v[i] = rad2deg * mod_angle_c(dlon + deg2rad * zone_number_to_central_longitude(zone))

    sample = easting if getattr(easting, "ndim", None) is not None else northing
    return (
        _maybe_numpy_array(sample, lat_arr),
        _maybe_numpy_array(sample, lon_arr),
    )

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
