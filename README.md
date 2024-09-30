[![test](https://img.shields.io/github/v/release/murphycj/agfusion)](https://img.shields.io/github/v/release/murphycj/agfusion)
[![test](https://img.shields.io/github/last-commit/murphycj/agfusion)](https://img.shields.io/github/last-commit/murphycj/agfusion)
[![Downloads](https://pepy.tech/badge/agfusion)](https://pepy.tech/project/agfusion)
[![Unit tests](https://github.com/murphycj/AGFusion/actions/workflows/test.yaml/badge.svg)](https://github.com/murphycj/AGFusion/actions/workflows/test.yaml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F080903-lightgrey.svg?style=flat-square)](http://dx.doi.org/10.1101/080903)
[![Citations: 21](https://img.shields.io/badge/citations-21-brightgreen)](https://scholar.google.com/scholar?cites=10329835253652722473&as_sdt=5,33&sciodt=0,33&hl=en)


# Annotate Gene Fusion (AGFusion)
**Checkout the webapp:** [https://www.agfusion.app](https://www.agfusion.app)

AGFusion (pronounced 'A G Fusion') is a python package for annotating gene fusions from the human or mouse genomes. AGFusion simply needs the reference genome, the two gene partners, and the fusion junction coordinates as input, and outputs the following:

* FASTA files of cDNA, CDS, and protein sequences.
* Visualizes the protein domain and exon architectures of the fusion transcripts.
* Saves tables listing the coordinates of protein features and exons included in the fusion.
* Optional exon structure and protein domain visualization of the wild-type  version of the fusion gene partners.

Some other things to know:

* AGFusion automatically predicts the functional effect of the gene fusion (e.g. in-frame, out-of-frame, etc.).
* Annotation is by default done only for canonical gene isoforms, but there is the option to annotate all gene non-canonical isoform combinations.
* All gene and protein annotation is from Ensembl
* Supports up to Ensembl release 95


## Table of Contents

- [Installation](#installation)
- [Dependencies](#dependencies)
- [Examples](#examples)
  * [Basic Usage](#basic-usage)
  * [Plotting wild-type protein and exon structure](#plotting-wild-type-protein-and-exon-structure)
  * [Canonical gene isoforms](#canonical-gene-isoforms)
  * [Input from fusion-finding algorithms](#input-from-fusion-finding-algorithms)
  * [Graphical parameters](#graphical-parameters)
- [Building your own database](#building-your-own-database)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Citing AGFusion](#citing-agfusion)

# Installation

**Step 1:** Install AGFusion.

```
pip install agfusion
```

**Step 2:** Download your desired pyensembl reference genome database. For example:

```
For GRCh38/hg38:
pyensembl install --species homo_sapiens --release 95

For GRCh37/hg19:
pyensembl install --species homo_sapiens --release 75

For GRCm38/mm10:
pyensembl install --species mus_musculus --release 87
```

**Step 3:** Finally, download your desired AGFusion database.

```
For GRCh38/hg38:
agfusion download -g hg38

For GRCh37/hg19:
agfusion download -g hg19

For GRCm38/mm10:
agfusion download -g mm10
```

You can view all supported species and ensembl releases with ```agfusion download -a```.

# Dependencies

* Python 3.7 or higher
* Python package dependencies are listed in [requirements.txt](https://github.com/murphycj/AGFusion/blob/master/requirements.txt).

## Examples

### Basic Usage

You just need to provide the two fusion gene partners (gene symbol, Ensembl ID, or Entrez gene ID), their predicted fusion junctions in genomic coordinates, and the genome build. You can also specify certain transcripts with Ensembl transcript ID or RefSeq ID

Example usage from the command line:

```
agfusion annotate \
  --gene5prime DLG1 \
  --gene3prime BRAF \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  -db agfusion.mus_musculus.87.db \
  -o DLG1-BRAF
```

The protein domain structure of the DLG1-BRAF fusion:

![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000064477-ENSMUST00000002487.png)

The exon structure of the DLG1-BRAF fusion:

![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000064477-ENSMUST00000002487.exon.png)

### Plotting wild-type protein and exon structure

You can additionally plot the wild-type proteins and exon structures for each gene with --WT flag.

```
agfusion annotate \
   -g5 ENSMUSG00000022770 \
   -g3 ENSMUSG00000002413 \
   -j5 31684294 \
   -j3 39648486 \
   -db agfusion.mus_musculus.87.db \
   -o DLG1-BRAF \
   --WT
```

### Canonical gene isoforms

By default AGFusion only plots the [canonical](http://useast.ensembl.org/Help/Glossary?id=346) gene isoforms, but you can tell AGFusion to include non-canonical isoform with the --noncanonical flag.

```
agfusion annotate \
  -g5 ENSMUSG00000022770 \
  -g3 ENSMUSG00000002413 \
  -j5 31684294 \
  -j3 39648486 \
  -db agfusion.mus_musculus.87.db \
  -o DLG1-BRAF \
  --noncanonical
```

### Input from fusion-finding algorithms

You can provide as input output files from fusion-finding algorithms. Currently supported algorithms are:

* [Arriba](https://github.com/suhrig/arriba)
* Bellerophontes
* BreakFusion
* ChimeraScan
* ChimeRScope
* [deFuse](https://github.com/amcpherson/defuse)
* EricScript
* [FusionCatcher](https://github.com/ndaniel/fusioncatcher)
* FusionHunter
* FusionMap
* InFusion
* [JAFFA](https://github.com/Oshlack/JAFFA)
* [LongGF](https://github.com/WGLab/LongGF)
* MapSplice (only if --gene-gtf specified)
* [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion)
* TopHat-Fusion



Below is an example for FusionCatcher.

```
agfusion batch \
  -f final-list_candidate-fusion-genes.txt \
  -a fusioncatcher \
  -o test \
  -db agfusion.mus_musculus.87.db
```

### Graphical parameters

You can change domain names and colors:

```
agfusion annotate \
  -g5 ENSMUSG00000022770 \
  -g3 ENSMUSG00000002413 \
  -j5 31684294 \
  -j3 39648486 \
  -db agfusion.mus_musculus.87.db \
  -o DLG1-BRAF \
  --recolor "Pkinase_Tyr;red" --recolor "L27_1;blue" \
  --rename "Pkinase_Tyr;Kinase" --rename "L27_1;L27"
```

![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000064477-ENSMUST00000002487-recolorRename.png)

You can rescale the protein length so that images of two different fusions have appropriate relative lengths when plotted side by side:

```
agfusion annotate \
  -g5 ENSMUSG00000022770 \
  -g3 ENSMUSG00000002413 \
  -j5 31684294 \
  -j3 39648486 \
  -db agfusion.mus_musculus.87.db \
  -o DLG1-BRAF \
  --recolor "Pkinase_Tyr;red" --recolor "L27_1;blue" \
  --rename "Pkinase_Tyr;Kinase" --rename "L27_1;L27" \
  --scale 2000
agfusion annotate \
  -g5 FGFR2 \
  -g3 DNM3 \
  -j5 130167703 \
  -j3 162019992 \
  -db agfusion.mus_musculus.87.db \
  -o FGFR2-DNM3 \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2000
```

![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000064477-ENSMUST00000002487-rescale.png)
![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000122054-ENSMUST00000070330-rescale.png)

# Building your own database
AGFusion uses a pre-built SQLite database to annotation gene fusions; in addition to data from pyensembl. The SQLite databases are stored on AWS S3.

Follow the steps below if you want to build your own SQLite database:

(1) Install [mysqlclient](https://github.com/PyMySQL/mysqlclient).

(2) Download and unzip the PFAM reference file: [https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz)

(3) Install your desired pyensembl reference genome. For example: `pyensembl install --release 111`.

(4) Build the AGFusion database: `agfusion build -d . -s homo_sapiens -r 111 --pfam Pfam-A.clans.tsv`

# Troubleshooting

**(1) Problem:** I get a warning message like the following:


> 2017-08-28 15:02:51,377 - AGFusion - WARNING - No cDNA sequence available for AC073283.4! Will not print cDNA sequence for the AC073283.4-MSH2 fusion. You might be working with an outdated pyensembl. Update the package and rerun 'pyensembl install'

**Solution:** Run the following to update pyensembl package and database:

```
git clone git@github.com:hammerlab/pyensembl.git
cd pyensembl
sudo pip install .
pyensembl install --release (your-release) --species (your-species)
```

**(2) Problem:** Cannot run `agfusion download` due to `URLError`.
When downloading the database you may run into this error:
```
urllib.error.URLError: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate (_ssl.c:1108)>
```

**Solution:**
A potential solution for Mac users is from [here](https://stackoverflow.com/questions/27835619/urllib-and-ssl-certificate-verify-failed-error). You can run the following command:

```
/Applications/Python\ 3.8/Install\ Certificates.command
```


# License

MIT license

# Citing AGFusion

You can cite bioRxiv: http://dx.doi.org/10.1101/080903
