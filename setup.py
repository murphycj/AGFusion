"""
AGFusion is a python package to annotate and visualize gene fusions.
"""

import re

from setuptools import find_packages, setup

VERSIONFILE = "agfusion/version.py"
version = open(VERSIONFILE, "rt", encoding="utf8").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, version, re.M)
if mo:
    version = mo.group(1)
else:
    raise RuntimeError(f"Unable to find version string in {VERSIONFILE}.")

README = open("README.md", encoding="utf8").read()

setup(
    # version=version,
    version="1.3.0",
    name="agfusion",
    packages=find_packages(),
    description="Python package to annotate and visualize gene fusions.",
    author="Charles Murphy",
    author_email="murphy.charlesj@gmail.com",
    license="MIT",
    url="https://github.com/murphycj/AGFusion",
    long_description=README,
    long_description_content_type="text/markdown",
    include_package_data=True,
    scripts=["bin/agfusion"],
    install_requires=[
        "matplotlib>=1.5.0",
        "pandas>=0.18.1",
        "biopython>=1.67",
        "future>=0.16.0",
        "pyensembl>=1.1.0",
    ],
)
