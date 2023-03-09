"""linz-geodetic setup module
"""

# Comments from samply python packaging script

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, "DESCRIPTION.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="linz-geodetic",
    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version="1.0.0",
    description="LINZ.geodetic module - generic geodetic functions",
    long_description=long_description,
    # The project's main homepage.
    url="https://github.com/linz/python-linz-geodetic",
    # Author details
    author="Chris Crook",
    author_email="ccrook@linz.govt.nz",
    # Choose your license
    license="MIT",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Stable",
        # Indicate who your project is intended for
        "Intended Audience :: LINZ coders",
        "Topic :: Geodetic",
        # Pick your license as you wish (should match "license" above)
        "License :: OSI Approved :: GPL License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        "Programming Language :: Python :: 3",
    ],
    # What does your project relate to?
    keywords="geodesy",
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),
    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=["numpy"],
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={},
    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={},
    # Namespace package - other modules may include into these packages
    namespace_packages=["LINZ", "LINZ.Geodetic"],
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        "console_scripts": [
            "ellipsoid=LINZ.Geodetic.Ellipsoid:main",
        ],
    },
)
