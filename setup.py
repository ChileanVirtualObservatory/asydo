import sys
from setuptools import setup

def get_dependencies():
    if sys.version_info.major == 2:
        return ['numpy', 'astropy', 'matplotlib', 'urllib3',
         'scipy>=0.19.0', 'pysqlite']
    else:
        return ['numpy', 'astropy', 'matplotlib', 'urllib3',
                'scipy>=0.19.0']

def get_python_classifier():
    if sys.version_info.major == 2 and sys.version_info.minor >= 7:
        return 'Programming Language :: Python :: 2.7'
    else:
        return 'Programming Language :: Python :: 3.6'

setup(
    name = 'asydo',
    version = '0.1.0',
    description = 'Astronomical Syntetic Data Observations',
    long_description='''
            We propose an Astronomical SYnthetic Data Observatory (ASYDO),
            a virtual service that generates synthetic spectroscopic data in
            the form of data cubes. The objective of the tool is not to produce
            accurate astrophysical simulations, but to generate a large number
            of labelled synthetic data, to assess advanced computing algorithms
            for astronomy and to develop novel Big Data algorithms.''',
    url = 'https://github.com/ChileanVirtualObservatory/asydo',
    author = 'CSRG',
    author_email = 'contact@lirae.cl',
    license = 'GPLv2',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Topic :: Scientific/Engineering :: Astronomy',
        get_python_classifier(),
    ],
    zip_safe = True,
    package_dir = {'': 'src'},
    packages = ['asydo'],
    install_requires = get_dependencies()
)
