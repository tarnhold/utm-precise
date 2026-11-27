# cython: language_level=3, boundscheck=False, wraparound=False
# cython: cdivision=True, initializedcheck=False

from libc.math cimport (
    sin, cos, tan, atan, atan2, sinh, cosh, asinh, atanh,
    sqrt, hypot, exp, fabs, fmax, fmod, M_PI,
)

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
# first eccentricity squared (e^2) and first eccentricity e
cdef double E2 = (2.0 - F) * F
cdef double ES = sqrt(E2)
cdef double E2M = 1.0 - E2
# third flattening n
cdef double N = F / (2.0 - F)

# Rectified-meridian scale A = a * b1, with
# b1 = (1 + n^2/4 + n^4/64 + n^6/256) / (1+n)
# (Karney 2011 / Poder-Engsager, 6th-order series in n)
cdef double _n2 = N * N
cdef double B1 = (_n2 * _n2 * _n2 + 4.0 * _n2 * _n2 + 64.0 * _n2 + 256.0) / (
    256.0 * (1.0 + N)
)
cdef double A = R * B1
cdef double AK0 = A * K0

# Forward Krueger coefficients alpha[1..6] (Karney 2011, order 6)
cdef double ALP1 = N * (
    ((((31564.0 * N - 66675.0) * N + 34440.0) * N + 47250.0) * N - 100800.0) * N
    + 75600.0
) / 151200.0
cdef double ALP2 = _n2 * (
    (((-1983433.0 * N + 863232.0) * N + 748608.0) * N - 1161216.0) * N + 524160.0
) / 1935360.0
cdef double ALP3 = _n2 * N * (
    ((670412.0 * N + 406647.0) * N - 533952.0) * N + 184464.0
) / 725760.0
cdef double ALP4 = _n2 * _n2 * (
    (6601661.0 * N - 7732800.0) * N + 2230245.0
) / 7257600.0
cdef double ALP5 = _n2 * _n2 * N * (-13675556.0 * N + 3438171.0) / 7983360.0
cdef double ALP6 = _n2 * _n2 * _n2 * 212378941.0 / 319334400.0

# Inverse Krueger coefficients beta[1..6]
cdef double BET1 = N * (
    ((((384796.0 * N - 382725.0) * N - 6720.0) * N + 932400.0) * N - 1612800.0) * N
    + 1209600.0
) / 2419200.0
cdef double BET2 = _n2 * (
    (((-1118711.0 * N + 1695744.0) * N - 1174656.0) * N + 258048.0) * N + 80640.0
) / 3870720.0
cdef double BET3 = _n2 * N * (
    ((22276.0 * N - 16929.0) * N - 15984.0) * N + 12852.0
) / 362880.0
cdef double BET4 = _n2 * _n2 * ((-830251.0 * N - 158400.0) * N + 197865.0) / 7257600.0
cdef double BET5 = _n2 * _n2 * N * (-435388.0 * N + 453717.0) / 15966720.0
cdef double BET6 = _n2 * _n2 * _n2 * 20648693.0 / 638668800.0

cdef double TAU_TOL = 1.4901161193847657e-9   # sqrt(eps) / 10
cdef double TAU_MAX = 134217728.0            # 2 / sqrt(eps)
cdef double ES_ATANH_ES = ES * atanh(ES)

ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"


cdef inline double mod_angle_c(double value) noexcept:
    """Angle in radians in range [-pi, pi)."""
    cdef double rem = fmod(value + M_PI, 2.0 * M_PI)
    if rem < 0.0:
        rem += 2.0 * M_PI
    return rem - M_PI


cpdef mod_angle(double value):
    """Returns angle in radians in range [-pi, pi)."""
    return mod_angle_c(value)


cdef inline double taupf(double tau) noexcept:
    """tan(chi) from tan(phi), chi = conformal latitude (Karney)."""
    cdef double tau1 = hypot(1.0, tau)
    cdef double sig = sinh(ES * atanh(ES * tau / tau1))
    return hypot(1.0, sig) * tau - sig * tau1


cdef inline double tauf(double taup) noexcept:
    """tan(phi) from tan(chi) by Newton's method (Karney)."""
    cdef double tau, taupa, dtau, stol
    cdef int i
    if fabs(taup) > 70.0:
        tau = taup * exp(ES_ATANH_ES)
    else:
        tau = taup / E2M
    if not (fabs(tau) < TAU_MAX):
        return tau
    stol = TAU_TOL * fmax(1.0, fabs(taup))
    for i in range(8):
        taupa = taupf(tau)
        dtau = (taup - taupa) * (1.0 + E2M * tau * tau) / (
            E2M * hypot(1.0, tau) * hypot(1.0, taupa)
        )
        tau += dtau
        if fabs(dtau) < stol:
            break
    return tau


cdef void tm_forward(double lat_rad, double dlon, double *x, double *y) noexcept:
    """Transverse Mercator forward (metres from CM / equator, no false offsets)."""
    cdef double tau, taup, c_lam, s_lam, xi_p, eta_p, xi, eta
    cdef double s2, s4, s6, s8, s10, s12
    cdef double c2, c4, c6, c8, c10, c12
    cdef double sh2, sh4, sh6, sh8, sh10, sh12
    cdef double ch2, ch4, ch6, ch8, ch10, ch12

    tau = tan(lat_rad)
    taup = taupf(tau)
    c_lam = cos(dlon)
    s_lam = sin(dlon)
    xi_p = atan2(taup, c_lam)
    eta_p = asinh(s_lam / sqrt(taup * taup + c_lam * c_lam))

    s2 = sin(2.0 * xi_p)
    c2 = cos(2.0 * xi_p)
    # Chebyshev recurrence: sin((k+1)th) = 2 cos(2ξ') sin(kth) - sin((k-1)th)
    s4 = 2.0 * c2 * s2
    c4 = 2.0 * c2 * c2 - 1.0
    s6 = 2.0 * c2 * s4 - s2
    c6 = 2.0 * c2 * c4 - c2
    s8 = 2.0 * c2 * s6 - s4
    c8 = 2.0 * c2 * c6 - c4
    s10 = 2.0 * c2 * s8 - s6
    c10 = 2.0 * c2 * c8 - c6
    s12 = 2.0 * c2 * s10 - s8
    c12 = 2.0 * c2 * c10 - c8

    sh2 = sinh(2.0 * eta_p)
    ch2 = cosh(2.0 * eta_p)
    sh4 = 2.0 * ch2 * sh2
    ch4 = 2.0 * ch2 * ch2 - 1.0
    sh6 = 2.0 * ch2 * sh4 - sh2
    ch6 = 2.0 * ch2 * ch4 - ch2
    sh8 = 2.0 * ch2 * sh6 - sh4
    ch8 = 2.0 * ch2 * ch6 - ch4
    sh10 = 2.0 * ch2 * sh8 - sh6
    ch10 = 2.0 * ch2 * ch8 - ch6
    sh12 = 2.0 * ch2 * sh10 - sh8
    ch12 = 2.0 * ch2 * ch10 - ch8

    xi = (xi_p
          + ALP1 * s2 * ch2
          + ALP2 * s4 * ch4
          + ALP3 * s6 * ch6
          + ALP4 * s8 * ch8
          + ALP5 * s10 * ch10
          + ALP6 * s12 * ch12)
    eta = (eta_p
           + ALP1 * c2 * sh2
           + ALP2 * c4 * sh4
           + ALP3 * c6 * sh6
           + ALP4 * c8 * sh8
           + ALP5 * c10 * sh10
           + ALP6 * c12 * sh12)
    x[0] = AK0 * eta
    y[0] = AK0 * xi


cdef void tm_reverse(double x, double y, double *lat_rad, double *dlon) noexcept:
    """Transverse Mercator inverse (metres from CM / equator, no false offsets)."""
    cdef double xi, eta, xi_p, eta_p, s, c, r, taup
    cdef double s2, s4, s6, s8, s10, s12
    cdef double c2, c4, c6, c8, c10, c12
    cdef double sh2, sh4, sh6, sh8, sh10, sh12
    cdef double ch2, ch4, ch6, ch8, ch10, ch12

    xi = y / AK0
    eta = x / AK0

    s2 = sin(2.0 * xi)
    c2 = cos(2.0 * xi)
    s4 = 2.0 * c2 * s2
    c4 = 2.0 * c2 * c2 - 1.0
    s6 = 2.0 * c2 * s4 - s2
    c6 = 2.0 * c2 * c4 - c2
    s8 = 2.0 * c2 * s6 - s4
    c8 = 2.0 * c2 * c6 - c4
    s10 = 2.0 * c2 * s8 - s6
    c10 = 2.0 * c2 * c8 - c6
    s12 = 2.0 * c2 * s10 - s8
    c12 = 2.0 * c2 * c10 - c8

    sh2 = sinh(2.0 * eta)
    ch2 = cosh(2.0 * eta)
    sh4 = 2.0 * ch2 * sh2
    ch4 = 2.0 * ch2 * ch2 - 1.0
    sh6 = 2.0 * ch2 * sh4 - sh2
    ch6 = 2.0 * ch2 * ch4 - ch2
    sh8 = 2.0 * ch2 * sh6 - sh4
    ch8 = 2.0 * ch2 * ch6 - ch4
    sh10 = 2.0 * ch2 * sh8 - sh6
    ch10 = 2.0 * ch2 * ch8 - ch6
    sh12 = 2.0 * ch2 * sh10 - sh8
    ch12 = 2.0 * ch2 * ch10 - ch8

    xi_p = (xi
            - BET1 * s2 * ch2
            - BET2 * s4 * ch4
            - BET3 * s6 * ch6
            - BET4 * s8 * ch8
            - BET5 * s10 * ch10
            - BET6 * s12 * ch12)
    eta_p = (eta
             - BET1 * c2 * sh2
             - BET2 * c4 * sh4
             - BET3 * c6 * sh6
             - BET4 * c8 * sh8
             - BET5 * c10 * sh10
             - BET6 * c12 * sh12)

    s = sinh(eta_p)
    c = cos(xi_p)
    r = hypot(s, c)
    dlon[0] = atan2(s, c)
    if r == 0.0:
        lat_rad[0] = M_PI / 2.0 if xi_p >= 0.0 else -M_PI / 2.0
        return
    taup = sin(xi_p) / r
    lat_rad[0] = atan(tauf(taup))


cdef inline int zone_cm_deg(int zone_number) noexcept:
    return (zone_number - 1) * 6 - 180 + 3


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
    return zone_cm_deg(zone_number)


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
    """Convert latitude and longitude (GRS80, degrees) to UTM.

    Parameters
    ----------
    latitude, longitude : float or array-like
        Geodetic coordinates in degrees. Longitude in [-180, 180).
        Array-like inputs (lists or NumPy arrays) are converted element-wise.
        NumPy is optional; arrays work without it and return lists, or
        ndarrays if a NumPy array was passed in.
    force_zone_number : int or array-like, optional
        Force conversion into this UTM zone.
    strict : bool
        If True (default), reject latitudes outside 80°S–84°N.
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
            raise OutOfRangeError(
                'latitude out of range (must be between 80 deg S and 84 deg N)'
            )
    if not -180.0 <= longitude < 180.0:
        raise OutOfRangeError(
            'longitude out of range (must be between 180 deg W and 180 deg E)'
        )

    if force_zone_number is None:
        zone_number = latlon_to_zone_number(latitude, longitude)
    else:
        zone_number = force_zone_number

    if not 1 <= zone_number <= 60:
        raise OutOfRangeError('zone number out of range (must be between 1 and 60)')

    zone_letter = latitude_to_zone_letter(latitude)
    dlon = mod_angle_c(deg2rad * (longitude - zone_cm_deg(zone_number)))
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
        dlon = mod_angle_c(deg2rad * (lon - zone_cm_deg(zone_number)))
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
    """Convert a UTM coordinate to latitude and longitude (GRS80, degrees).

    Parameters
    ----------
    easting, northing : float or array-like
        UTM coordinates in metres. Array-like inputs are converted
        element-wise. NumPy is optional.
    zone_number : int or array-like
        UTM zone number (1–60).
    zone_letter : str, optional
        UTM zone letter. Either this or ``northern`` must be given.
    northern : bool, optional
        True for the northern hemisphere, False for the southern.
    strict : bool
        If True (default), reject easting/northing outside the usual UTM range.
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
    cdef double x, y, lat, dlon, lon0
    cdef bint north

    north = _hemisphere(zone_letter, northern, zone_number)
    if strict:
        if not 100000.0 <= easting < 1000000.0:
            raise OutOfRangeError(
                'easting out of range (must be between 100.000 m and 999.999 m)'
            )
        if not 0.0 <= northing <= 10000000.0:
            raise OutOfRangeError(
                'northing out of range (must be between 0 m and 10.000.000 m)'
            )
    if not 1 <= zone_number <= 60:
        raise OutOfRangeError('zone number out of range (must be between 1 and 60)')

    x = easting - 500000.0
    y = northing
    if not north:
        y -= 10000000.0
    tm_reverse(x, y, &lat, &dlon)
    lon0 = deg2rad * zone_cm_deg(zone_number)
    return rad2deg * lat, rad2deg * mod_angle_c(dlon + lon0)


cdef bint _hemisphere(object zone_letter, object northern, int zone_number):
    if zone_letter:
        zone_letter = zone_letter.upper()
        if len(zone_letter) > 1:
            raise OutOfRangeError('zone letter out of range (must be one character)')
        if not 'C' <= zone_letter <= 'X' or zone_letter in ['I', 'O']:
            raise OutOfRangeError('zone letter out of range (must be between C and X)')
        if zone_letter == 'X' and zone_number in [32, 34, 36]:
            raise OutOfRangeError(
                'zone number not valid (there is no zone ' + str(zone_number) + ' in X)'
            )
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
        lon_v[i] = rad2deg * mod_angle_c(dlon + deg2rad * zone_cm_deg(zone))

    sample = easting if getattr(easting, "ndim", None) is not None else northing
    return (
        _maybe_numpy_array(sample, lat_arr),
        _maybe_numpy_array(sample, lon_arr),
    )
