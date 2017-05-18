# Annotate Gene Fusion (AGFusion)
For a given gene fusion, AGFusion will predict the cDNA, CDS, and protein sequences resulting from fusion of all combinations of transcripts and save them to fasta files. AGFusion can also plot the protein domain architecture of the fusion transcripts. Currently, only PFAM domains are used to annotate gene fusions. CDS and protein sequences are only outputted for fusions that can form a protein.

Docs are at http://pythonhosted.org/agfusion/

## Table of Contents

- [Examples](#Examples)
  * [Basic Usage](#Basic-Usage)
  * [Plotting wild-type protein and exon structure](#Plotting-wild-type-protein)
  * [Canonical gene isoforms](#Canonical-gene-isoforms)
  * [Input from fusion-finding algorithms](#Input-from-fusion-finding-algorithms)
  * [Domain names and colors](#Domain-names-and-colors)
  * [Re-scaling protein length](#Re-scaling-protein-length)
- [Installation](#Installation)
- [Dependencies](#Dependencies)
- [License](#License)
- [Citing AGFusion](#Citing-AGFusion)


## Examples

### Basic Usage

You just need to provide the two fusion gene partners (gene symbol, Ensembl ID, or Entrez gene ID), their predicted fusion junctions in genomic coordinates, and the genome build. You can also specify certain transcripts with Ensembl transcript ID or RefSeq ID

Example usage from the command line:

```
agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF
```

The protein domain structure of the DLG1-BRAF fusion:

![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000064477-ENSMUST00000002487.png)

The exon structure of the DLG1-BRAF fusion:

![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000064477-ENSMUST00000002487.exon.png)

The cDNA, CDS, and protein fusion sequences are outputted in fasta format. Here is an example protein sequence for the in-frame fusion:

```
>ENSMUSP00000064280-ENSMUSP00000002487 length=535, kD: 60.0739853, transcripts: ENSMUST00000064477-ENSMUST00000002487, genes: Dlg1-Braf, effect: in-frame
MPVRKQDTQRALHLLEEYRSKLSQTEDRQLRSSIERVINIFQSNLFQALIDIQEFYEVTL
LDNPKCVDHSKQCEPVQPVTTWEIASLPSTAVTSETLPGSLSPPVEDLIRDQGFRGDGAP
LNQLMRCLRKYQSRTPSPLLHSVPSEIVFDFEPGPVFRGSTTGLSATPPASLPGSLTNVK
ALQKSPGPQRERKSSSSSSSEDRSRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVY
KGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWC
EGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKI
GDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTG
QLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILA
SIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGEFAAFK
```

### Plotting wild-type protein and exon structure

You can additionally plot the wild-type proteins for each gene with --WT flag.

```
agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF \
  --WT
```

### Canonical gene isoforms

By default AGFusion only plots the [canonical](http://useast.ensembl.org/Help/Glossary?id=346) gene isoforms, but you can tell AGFusion to include non-canonical isoform with the --noncanonical flag. This can produce many fusion isoform combinations.

```
agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF \
  --noncanonical
```

### Input from fusion-finding algorithms

You can provide as input output files from fusion-finding algorithms. Currently available algorithms are: FusionCatcher, TopHat-Fusion, and STAR-Fusion.
Below is an example for FusionCatcher.

```
agfusion batch \
  --file final-list_candidate-fusion-genes.txt \
  -a fusioncatcher \
  -o test \
  -g GRCm38
```

### Domain names and colors

You can change domain names and colors:

```
agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF \
  --recolor "Pkinase_Tyr;red" --recolor "L27_1;blue" \
  --rename "Pkinase_Tyr;Kinase" --rename "L27_1;L27"
```

![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000064477-ENSMUST00000002487-recolorRename.png)

### Re-scaling protein length

You can rescale the protein length so that images of two different fusions have appropriate relative lengths when plotted side by side:

```
agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF \
  --recolor "Pkinase_Tyr;red" --recolor "L27_1;blue" \
  --rename "Pkinase_Tyr;Kinase" --rename "L27_1;L27" \
  --scale 2000
agfusion annotate \
  --gene5prime FGFR2 \
  --gene3prime DNM3 \
  --junction5prime 130167703 \
  --junction3prime 162019992 \
  --genome GRCm38 \
  --out FGFR2-DNM3 \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2000
```

![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000064477-ENSMUST00000002487-rescale.png)
![alt tag](https://github.com/murphycj/AGFusion/blob/master/doc/ENSMUST00000122054-ENSMUST00000070330-rescale.png)

# Installation

First you need to install pyensembl (and the other dependencies listed at the bottom) and download the reference genome you will use by running one of the following.

For GRCh38/hg38:

```
pyensembl install --release 84 --species homo_sapiens
```

For GRCh37/hg19:

```
pyensembl install --release 75 --species homo_sapiens
```

For GRCm38/mm10:

```
pyensembl install --release 84 --species mus_musculus
```

Then you can install AGFusion via the following:

```
pip install agfusion
```

Finally, download the AGFusion database:

```
agfusion download
```

# Dependencies

- python 2.7, 3.5
- matplotlib>=1.5.0
- pandas>=0.18.1
- biopython>=1.67
- future>=0.16.0
- pyensembl>=1.1.0

# License

MIT license

# Citing AGFusion

Manuscript under review. You can cite bioRxiv for now: http://dx.doi.org/10.1101/080903
