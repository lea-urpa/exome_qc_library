from setuptools import setup, find_packages

VERSION = '2.0.0'
DESCRIPTION = 'Exome seq QC with Hail'
LONG_DESCRIPTION = 'Exome sequencing quality control pipeline with Hail'

# Run with python setup.py sdist bdist_wheel to crete source distribution and wheel for PyPi

setup(
    # Must match name of module
    name="exome_qc",
    version=VERSION,
    author="Lea M. Urpa",
    author_email="lea.urpa@helsinki.fi",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=[],

    keywords=['genetics', 'genomics', 'hail', 'sequencing', 'exome'],
    classifiers= [
        "Development Status :: 3 - Alpha", #TODO check the definitions on this
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3"
    ]
)