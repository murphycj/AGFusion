.. AGFusion documentation master file, created by
   sphinx-quickstart on Mon Oct 10 10:34:12 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AGFusion's documentation!
====================================

:Author: Charlie Murphy
:Contact: murphy.charlesj@gmail.com
:Source code: `GitHub <https://github.com/murphycj/AGFusion>`_
:License: MIT
:Packages: `PyPI <https://pypi.python.org/pypi/agfusion>`_
:Article: Pending...

For a given gene fusion, AGFusion will predict the cDNA, CDS, and protein sequences
resulting from fusion of all combinations of transcripts and save them to fasta
files. AGFusion can also plot the protein domain architecture of the fusion
transcripts. Currently, only PFAM domains are used to annotate gene fusions.
CDS and protein sequences are only outputted for fusions that can form a protein.

.. toctree::

    quickstart

.. gallery?

Command line usage
------------------

.. toctree::
    :maxdepth: 2

    examples

Python API
----------

.. toctree::
    :maxdepth: 2

    agfusion


Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
