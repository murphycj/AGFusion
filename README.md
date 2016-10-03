# Visualize Gene Fusion (VGFusion)
Python package that visualizes and annotates gene fusions. Currently, the software can only visualize the PFAM domains of in-frame gene fusions. The software can also output fasta files of the predicted cDNA, CDS, and protein sequences resulting from fusion of all combinations of transcripts.

# Example Usage

The minimum amount of information you need to provide are the two fusion gene partners (gene symbol or Ensembl ID) and their respective predicted fusion junctions in genomic coordinates.

Example usage within a python script:

```
python
import agfusion

data = pyensembl.EnsemblRelease(84,'mouse')
db = agfusion.AGFusionDB('../data/agfusion.db')

dlg1 = agfusion.Gene(
    gene="ENSMUSG00000022770",
    junction=31684294,
    db=db,
    pyensembl_data=data
)

braf = agfusion.Gene(
    gene="ENSMUSG00000002413",
    junction=39648486,
    db=db,
    pyensembl_data=data
)

fusion = agfusion.model.Fusion(dlg1,braf,db=db,middlestar=False)

fusion.save_transcript_cdna('DLG1-BRAF_mouse')
fusion.save_transcript_cds('DLG1-BRAF_mouse')
fusion.save_proteins('DLG1-BRAF_mouse')
fusion.save_images('DLG1-BRAF_mouse')
```

Example usage from the command line:

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

Example output visualization of the domain structure of the DLG1-BRAF fusion:

![alt tag](https://github.com/murphycj/AGFusion/tree/master/test/DLG1-BRAF/ENSMUST00000023454-ENSMUST00000002487.png)

# Dependencies

python 2.7.8
pandas
biopython
pyensembl
matplotlib
mpld3
json

# Installation

```
pip install agfusion
```

After install pyensembl you need to install the reference genome you will use.

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

# License

MIT license