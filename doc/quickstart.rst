Quick start
===========

To run AGFusion on your own machine, keep reading.


Install AGFusion
----------------

First you need to install `pyensembl <https://github.com/hammerlab/pyensembl>`_
(and the other dependencies listed at the bottom of this readme) and download
the reference genome you will use by running one of the following.

For GRCh38::

    pyensembl install --release 84 --species homo_sapiens

For GRCh37::

    pyensembl install --release 75 --species homo_sapiens

For GRCm38::

    pyensembl install --release 84 --species mus_musculus

Those are the only reference genomes supported right now (but more will be
added later).

Then run::

    pip install agfusion

Basic usage
-----------

For basic usage of AGFusion from the command line you need to provide the
5' and 3' gene fusion partners (gene symbol or Ensembl gene ID), their
respective fusion junction coordinates, the reference genome, and the
output directory. Here is a simple example of a DLG1-BRAF fusion::

    agfusion \
      --gene5prime DLG1 \
      --gene3prime ENSMUSG00000002413 \
      --junction5prime 31684294 \
      --junction3prime 39648486 \
      --genome GRCm38 \
      --out DLG1-BRAF

.. image:: ENSMUST00000132176-ENSMUST00000002487.png
    :align: center
