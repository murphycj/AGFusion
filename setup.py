import os
import re
import shutil
import site

from setuptools import setup, find_packages

VERSIONFILE = "agfusion/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

README = open('README.md').read()

setup(
    version=verstr,
    name='agfusion',
    packages=find_packages(),
    description="Python package to annotate and visualize gene fusions.",
    author='Charles Murphy',
    author_email='murphy.charlesj@gmail.com',
    license='MIT',
    url='https://github.com/murphycj/AGFusion',
    long_description=README,
    long_description_content_type='text/markdown',
    include_package_data=True,
    scripts=['bin/agfusion'],
    install_requires=[
        'matplotlib>=1.5.0',
        'pandas>=0.18.1',
        'biopython>=1.67',
        'future>=0.16.0',
        'pyensembl>=1.1.0'
    ]
)
