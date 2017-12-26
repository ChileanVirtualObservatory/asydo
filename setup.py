from setuptools import setup, find_packages

setup(
    name = 'asydo',
    version = '0.1.0.dev1',
    description = 'Astronomical Syntetic Data Observations',
    url = 'https://github.com/ChileanVirtualObservatory/asydo',
    author = 'CSRG',
    author_email = 'contact@lirae.cl',
    classifiers = [
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6'
    ],
    zip_safe = True,
    package_dir = {'': 'src'},
    packages = ['asydo'],
    install_requires = ['acalib', 'numpy', 'astropy', 'matplotlib', 'urllib3',
                        'scipy>=0.19.0']
)
