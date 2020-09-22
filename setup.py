from setuptools import setup, find_packages
import pathlib

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="vortexfitting",
    version="1.0.1",
    description="A tool to locate and characterize vortices",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/guilindner/VortexFitting",
    author="Guilherme Anrain Lindner",
    keywords='vortex cfd fluid mechanics',
    author_email="lindner.guilherme@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        'Programming Language :: Python :: 3.8',
    ],
    
    packages= ['vortexfitting'],
    #package_dir={'vortexfitting': 'vortexfitting'},
    #package_data={'vortexfitting': ['data/*.*']},
    #include_package_data=True,
    install_requires=['netcdf4','matplotlib', 'numpy', 'scipy'],
    project_urls={  
        'Source': 'https://github.com/guilindner/VortexFitting',
    },
    entry_points={
        "console_scripts": [
            "vortexfitting = vortexfitting.__main__:main",
        ]
    },
)
