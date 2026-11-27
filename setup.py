from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        "utm/conversion.pyx",
        compiler_directives={"language_level": 3},
    ),
    scripts=["scripts/utm-converter"],
)
