# Annotate Gene Fusion (AGFusion)
For a given gene fusion, AGFusion will predict the cDNA, CDS, and protein sequences resulting from fusion of all combinations of transcripts and save them to fasta files. AGFusion can also plot the protein domain architecture of the fusion transcripts. Currently, only PFAM domains are used to annotate gene fusions. CDS and protein sequences are only outputted for fusions that can form a protein.

# Example Usage

You just need to provide the two fusion gene partners (gene symbol or Ensembl ID), their predicted fusion junctions in genomic coordinates, and the genome build.

Example usage from the command line:

```
./agfusion \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF
```

Example output visualization of the domain structure of the DLG1-BRAF fusion:

![alt tag](https://github.com/murphycj/AGFusion/blob/master/agfusion/data/ENSMUST00000132176-ENSMUST00000002487.png)

You can programmatically change domain names and colors:

```
../bin/agfusion \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF \
  --colors Pkinase_Tyr:red L27_1:#00cc00 \
  --rename Pkinase_Tyr:Kinase L27_1:L27
```

![alt tag](https://github.com/murphycj/AGFusion/blob/master/agfusion/data/ENSMUST00000132176-ENSMUST00000002487-color.png)

You can rescale the protein length so that images of two different fusions have appropriate relative lengths when plotted side by side:

```
./agfusion \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF \
  --scale 2000

./agfusion \
  --gene5prime FGFR2 \
  --gene3prime DNM3 \
  --junction5prime 130167703 \
  --junction3prime 162019992 \
  --genome GRCm38 \
  --out FGFR2-DNM3 \
  --scale 2000

```

![alt tag](https://github.com/murphycj/AGFusion/blob/master/agfusion/data/ENSMUST00000132176-ENSMUST00000002487-scale.png)
![alt tag](https://github.com/murphycj/AGFusion/blob/master/agfusion/data/ENSMUST00000120187-ENSMUST00000086074.png)

# Installation

First you need to install pyensembl (and the other dependencies) and download the reference genome you will use by running one of the following.

For GRCh38:

```
pyensembl install --release 84 --species homo_sapiens
```

For GRCh37:

```
pyensembl install --release 75 --species homo_sapiens
```

For GRCm38:

```
pyensembl install --release 84 --species mus_musculus
```

Then you can install AGFusion via the following:

```
pip install agfusion
```

# More examples

Example usage within a python script:

```
import agfusion
import pyensembl

data = pyensembl.EnsemblRelease(84,'mouse')
db = agfusion.AGFusionDB(‘/path/to/agfusion.db')


#construct the gene fusion

fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000022770",
        gene5primejunction=31684294,
        gene3prime="ENSMUSG00000002413",
        gene3primejunction=39648486,
        db=db,
        pyensembl_data=data
    )

#save the predicted cDNA, CDS, and protein sequences

fusion.save_transcript_cdna('DLG1-BRAF_mouse')
fusion.save_transcript_cds('DLG1-BRAF_mouse')
fusion.save_proteins('DLG1-BRAF_mouse')

#save the visualizations of the fusion domain structure

fusion.save_images('DLG1-BRAF_mouse')
```

For annotating with PFAM domains, AGFusion depends on a simple sqlite3 database (agfusion.db in the above examples). You can either download the database that comes with this package (recommended) or rebuild the database with the following command:

```
./build_db \
  --db agfusion.db \
  --genome GRCm38
```

And then specify the database with the --db flag:

```
./agfusion \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --db agfusion.db \
  --genome GRCm38 \
  --out DLG1-BRAF
```

# Dependencies

- python 2.7.8
- pyensembl>=0.9.5
- matplotlib>=1.5.0
- biomart>=0.9.0
- pandas>=0.18.1
- biopython>=1.67
- mpld3>=0.2
- jsonpickle>=0.9.

# License

MIT license