# Visualize Gene Fusion (VGFusion)
Python package providing that can visualize different annotations of a gene fusion.

# Example Usage

Within a python script:

‘’’python
import agfusion

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

# will return ['HLA-A']
gene_names = data.gene_names_at_locus(contig=6, position=29945884)

# get all exons associated with HLA-A
exon_ids  = data.exon_ids_of_gene_name('HLA-A')
‘’’

From the command line:

‘’’
python ./agfusion

‘’’

# Dependencies

python 2.7.8
pandas
biopython
pyensembl
matplotlib
mpld3
json

# Installation

‘’’
pip install agfusion
‘’’

After install pyensembl you need to install the reference genome you will use.

For GRCh38:

‘’’
pyensembl install --release 84 --species homo_sapiens
‘’’

For GRCh37:

‘’’
pyensembl install --release 75 --species homo_sapiens
‘’’

For GRCm38:

‘’’
pyensembl install --release 84 --species mus_musculus
‘’’

# License

MIT license