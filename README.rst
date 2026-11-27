utm-precise
===========

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://github.com/tarnhold/utm-precise/blob/master/LICENSE
.. image:: https://github.com/tarnhold/utm-precise/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/tarnhold/utm-precise/actions/workflows/ci.yml

Geodetic UTM/GRS80 projection for Python.

This is a modified version of the UTM/WGS84 converter
`Turbo87/utm <https://github.com/Turbo87/utm>`_. The original package uses a
lower-order series expansion (metre-level). ``utm-precise`` extends those
series for geodetic precision comparable in typical UTM zones to libraries
such as pyproj or GeographicLib, projects onto the **GRS80** ellipsoid
(the basis for ETRS89), and compiles the conversion routines with Cython.

The import name remains ``utm``. The distribution name on a package index
would be ``utm-precise``, so it does not collide with the upstream PyPI
package `utm <https://pypi.org/project/utm/>`_ (Turbo87/utm).

Modifications
-------------

* Expanded power series from metre to sub-millimetre precise projection
* Projection onto the GRS80 ellipsoid, as this is the basis for ETRS89 coordinates
* Cythonization for performance
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
    (294408.6629413165, 5628897.512984848, 32, 'U')

Convert UTM coordinates back to latitude and longitude. Either a zone letter
or ``northern=True/False`` must be given::

    >>> utm.to_latlon(294408.6629413165, 5628897.512984848, 32, 'U')
    (50.77534556000039, 6.083886670003538)
    >>> utm.to_latlon(294408.6629413165, 5628897.512984848, 32, northern=True)
    (50.77534556000039, 6.083886670003538)

A zone can be forced with ``force_zone_number`` (useful near zone boundaries).
Set ``strict=False`` to allow values slightly outside the official UTM range.

Command-line converter::

    $ utm-converter latlon 50.77534556 6.08388667
    294408.662941,5628897.512985,32,U

    $ utm-converter utm 294408.6629413165 5628897.512984848 32 U
    50.7753455598,6.0838866700

Accuracy
--------

Formulas follow Snyder (USGS Professional Paper 1395), Kelly (1986),
Hofmann-Wellenhof et al. (1994), and meridian-distance series after Deakin
and Weintritt — not Karney's or Krüger's complex series. Independent checks
against GeographicLib reference values shipped in the test suite and against
pyproj 3.8.0 / PROJ 9.8.1 (GRS80, November 2026) showed:

* Forward conversion vs. the GeographicLib test values: at most about
  **0.5 µm**.
* Forward conversion vs. pyproj in a standard 6°-wide UTM zone: typically
  nanometres to a few micrometres; about **3 µm** at the equator near the
  zone edge (~3° from the central meridian) in the sampled points.
* Wider meridians (Norway/Svalbard exceptions, or ``force_zone_number``):
  still sub-millimetre, but larger. Example: 60°N, 3°E in zone 32
  (~6° off the central meridian) differed by about **18 µm** forward and
  about **0.35 mm** inversely; a forced zone 32 at the equator and 3°E
  differed by about **0.43 mm** forward.
* The GRS80 vs. WGS84 ellipsoid choice is visible at this level (about
  **0.12 mm** at Aachen).

The "sub-millimetre" claim therefore holds comfortably inside a normal UTM
zone and only barely at large longitude differences from the central meridian.

Tests
-----

After installing the package (the Cython extension must be built)::

    python -m unittest discover -s test -v

License
-------

MIT. See ``LICENSE``.
