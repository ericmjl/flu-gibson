from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='FluGibson',
    packages=['FluGibson'],  # this must be the same as the name above
    version='1.3',
    description='A Python package for designing influenza reverse genetics \
                 primers using the seamless cloning methods (e.g. Gibson \
                 assembly, CPEC assembly).',
    long_description=long_description,
    author='Eric J. Ma',
    author_email='ericmajinglong@gmail.com',
    license='MIT',

    # use the URL to the github repo
    url='https://github.com/ericmjl/flu-gibson',
    keywords=['biology', 'molecular biology', 'cloning'],  # arbitrary keywords
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    requires=['networkx', 'pandas', 'matplotlib', 'biopython'],
    package_dir={'FluGibson': 'FluGibson'},
    package_data={'readme': ['README.md'],
                  'FluGibson': ['plasmid_backbones/*.fasta']},

)
