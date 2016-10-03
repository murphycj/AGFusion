import os
from setuptools import setup, find_packages

README = open('README.md').read()

setup(
    version = 0.1,
    name='agfusion',
    packages = find_packages(),
    description = "Python package providing that can visualize different annotations of a gene fusion.",
    author='Charles Murphy',
    license='MIT',
    long_description=README,
    install_requires=[
        'pyensembl>=0.9.5',
        'matplotlib>=1.5.0',
        'biomart>=0.9.0',
        'pandas>=0.18.1',
        'biopython>=1.67',
        'mpld3>=0.2',
        'jsonpickle>=0.9.3'
    ]
)
