"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    # $ pip install vortexfitting
    name='vortexfitting',
    version='0.1',
    description='A tool to locate and characterize vortices',  # Optional
    long_description=long_description,
    #url='https://github.com/pypa/sampleproject',  # Optional
    author='Guilherme Anrain Lindner',
    #author_email='@gmail.com',  # Optional

    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords='vortex cfd fluid mechanics',  # Optional
    package_dir={'': 'src'},
    packages=find_packages(where='src'),  # Required
    python_requires='>=3.6',

    install_requires=['netcdf4','matplotlib', 'numpy', 'scipy'],  # Optional

    #package_data={  # Optional
    #    'sample': ['package_data.dat'],
    #},

    project_urls={  # Optional
        'Source': 'https://github.com/guilindner/VortexFitting',
    },
)
