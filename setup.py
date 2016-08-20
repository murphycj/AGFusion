import os
from setuptools import setup, find_packages

README = open('README.md').read()

setup(
    version = 0.1,
    name='agfusion',
    packages = find_packages(),
    description = "Python package providing that can visualize different annotations of a gene fusion.",
    author='Charles Murphy',
    license='Non-commercial',
    long_description=README
)
