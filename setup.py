# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
exec(open("methplotlib/version.py").read())

setup(
    name="methplotlib",
    version=__version__,  # noqa: F821
    description="Plot nucleotide modification data from nanopore sequencing",
    long_description=open(path.join(here, "README.md")).read(),
    long_description_content_type="text/markdown",
    scripts=[
        "scripts/differential_modification",
        "scripts/allele_specific_modification",
    ],
    url="https://github.com/wdecoster/methplotlib",
    author="Wouter De Coster",
    author_email="decosterwouter@gmail.com",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    keywords="nanopore methylation modification",
    packages=find_packages(include=["methplotlib", "methplotlib.*"]),
    python_requires=">=3.9",
    install_requires=[
        "plotly>=4.9.0",
        "numpy>=1.16.5",
        "pandas>=0.23.4",
        "pyranges>=0.0.77",
        "fisher>=0.1.9",
        "scikit-learn",
        "pyfaidx",
        "biopython",
        "pysam",
    ],
    package_data={"methplotlib": []},
    package_dir={"methplotlib": "methplotlib"},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "methplotlib=methplotlib.methplotlib:main",
        ],
    },
)
