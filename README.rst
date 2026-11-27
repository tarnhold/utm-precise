utm-precise
===========

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://github.com/tarnhold/utm-precise/blob/master/LICENSE
.. image:: https://github.com/tarnhold/utm-precise/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/tarnhold/utm-precise/actions/workflows/ci.yml

Geodetic UTM/GRS80 projection for Python.

This is a modified version of the UTM/WGS84 converter
`Turbo87/utm <https://github.com/Turbo87/utm>`_. The original package uses a lower-order series expansion (metre-level).
``utm-precise`` uses Krüger's 6th-order series in the third flattening
(Karney 2011; Poder-Engsager) for nanometre-level accuracy, projects onto
the **GRS80** ellipsoid (the basis for ETRS89), and compiles the conversion
routines with Cython.

The import name remains ``utm``. The distribution name on a package index
would be ``utm-precise``, so it does not collide with the upstream PyPI
package `utm <https://pypi.org/project/utm/>`_ (Turbo87/utm).

Modifications
-------------

* Krüger series in *n* through order *n*:sup:`6` (Karney 2011; Poder-Engsager)
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
    (50.77534556, 6.083886670000001)
    >>> utm.to_latlon(294408.66294138716, 5628897.512984831, 32, northern=True)
    (50.77534556, 6.083886670000001)

A zone can be forced with ``force_zone_number`` (useful near zone boundaries).
Set ``strict=False`` to allow values slightly outside the official UTM range.

Command-line converter::

    $ utm-converter latlon 50.77534556 6.08388667
    294408.662941,5628897.512985,32,U

    $ utm-converter utm 294408.66294138716 5628897.512984831 32 U
    50.7753455598,6.0838866700

Accuracy
--------

The projection uses Krüger's series in the third flattening *n* through
order *n*\ :sup:`6` (Karney 2011; Poder-Engsager). Independent checks
against pyproj 3.8.0 / PROJ 9.8.1 (GRS80, November 2026) showed maximum
errors of a few nanometres, both in a normal UTM zone (≤3° from the
central meridian) and at 6–9° (Norway/Svalbard exceptions or
``force_zone_number``), well below 1 µm. The GRS80 vs. WGS84 ellipsoid
choice is much larger at this level (about 0.12 mm at Aachen).

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
