utm-precise
===========

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://github.com/tarnhold/utm-precise/blob/master/LICENSE
.. image:: https://github.com/tarnhold/utm-precise/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/tarnhold/utm-precise/actions/workflows/ci.yml

Geodetic UTM/GRS80 projection for Python.

This is a modified version of the UTM/WGS84 converter
`Turbo87/utm <https://github.com/Turbo87/utm>`_. The original package uses a lower-order series expansion (metre-level).
``utm-precise`` keeps those same formulae — meridian arc after Deakin/Weintritt,
footpoint latitude after Deakin, and Snyder/Kelly/Hofmann-Wellenhof power
series for easting, northing and the inverse — and extends the series to
higher order. It projects onto the **GRS80** ellipsoid (the basis for ETRS89)
and compiles the conversion routines with Cython.

The import name remains ``utm``. The distribution name on a package index
would be ``utm-precise``, so it does not collide with the upstream PyPI
package `utm <https://pypi.org/project/utm/>`_ (Turbo87/utm).

Modifications
-------------

* Higher-order power series: meridian arc through :math:`e^{14}` (Deakin 2006;
  Weintritt 2013; further coefficients by sympy), footpoint latitude through
  :math:`n^7` (Deakin 2012; further coefficients by sympy), and easting,
  northing and inverse through :math:`a^9`/:math:`a^{10}`, including the
  :math:`O(e'^2)` terms at :math:`a^7`/:math:`a^8` that the usual references
  leave out (coefficients by sympy, Taylor expansion of the complex meridian
  mapping)
* Projection onto the GRS80 ellipsoid, as this is the basis for ETRS89 coordinates
* Cythonization for performance
* Optional list/NumPy array input
* Extended testing (see below)
* Some minor bugfixes

Installation
------------

The package needs a C compiler and Cython (pulled in as a build dependency)::

    pip install .

Usage
-----

Convert geodetic coordinates (degrees, GRS80) to UTM::

    >>> import utm
    >>> utm.from_latlon(50.77534556, 6.08388667)
    (294408.66294138716, 5628897.512984831, 32, 'U')

Convert UTM coordinates back to latitude and longitude. Either a zone letter
or ``northern=True/False`` must be given::

    >>> utm.to_latlon(294408.66294138716, 5628897.512984831, 32, 'U')
    (50.775345560000005, 6.083886669999976)
    >>> utm.to_latlon(294408.66294138716, 5628897.512984831, 32, northern=True)
    (50.775345560000005, 6.083886669999976)

A zone can be forced with ``force_zone_number`` (useful near zone boundaries).
Set ``strict=False`` to allow values slightly outside the official UTM range.

Command-line converter::

    $ utm-converter latlon 50.77534556 6.08388667
    294408.662941,5628897.512985,32,U

    $ utm-converter utm 294408.66294138716 5628897.512984831 32 U
    50.7753455598,6.0838866700

Accuracy
--------

Independent checks against pyproj 3.8.0 / PROJ 9.8.1 (GRS80, November 2026).
The "original" column is the series as on ``master`` (meridian to :math:`e^{10}`,
footpoint to :math:`n^5`, easting/northing to :math:`a^7`/:math:`a^8` with
incomplete higher terms). "Extended" is this package. Krüger's 6th-order
series in the third flattening is shown only as a reference; this package
does not use it.

========================  ===========  ===========  ===========
Max. error                original     extended     Krüger n^6
========================  ===========  ===========  ===========
normal zone, forward      3.16 µm      3.3 nm       3.8 nm
normal zone, inverse      0.91 µm      3.2 nm       3.4 nm
6° from CM, forward       0.435 mm     0.15 µm      2.2 nm
6° from CM, inverse       0.353 mm     3.0 µm       2.0 nm
9° from CM, forward       8.27 mm      12.8 µm      2.2 nm
9° from CM, inverse       14.2 mm      0.27 mm      2.0 nm
========================  ===========  ===========  ===========

In a normal UTM zone (at most 3° from the central meridian) the extended
series is at the nanometre level. At 6–9° (Norway/Svalbard exceptions or
``force_zone_number``) the remaining truncation is far smaller than with
the original series, but the inverse at 9° is still in the tenths of a
millimetre. The GRS80 vs. WGS84 ellipsoid choice is much larger than the
normal-zone figure (about 0.12 mm at Aachen).

Array input
-----------

``from_latlon`` and ``to_latlon`` accept lists or, if installed, NumPy
arrays, and return lists or ``ndarray``\ s to match the input. NumPy is
an optional extra (``pip install utm-precise[numpy]``) and is not required
to use the library.

Tests
-----

After installing the package (the Cython extension must be built)::

    python -m unittest discover -s test -v

License
-------

MIT. See ``LICENSE``.
