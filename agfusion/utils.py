PROTEIN_DOMAIN = [
    ['prosite','prosite_start','prosite_end'],
    ['pfam','pfam_start','pfam_end'],
    ['smart','smart_start','smart_end'],
    ['pirsf','pirsf_start','pirsf_end'],
    ['superfamily','superfamily_start','superfamily_end'],
    ['hamap','hamap_start','hamap_end'],
    ['profile','profile_start','profile_end'],
    ['prints','prints_start','prints_end'],
    ['tigrfam','tigrfam_start','tigrfam_end'],
    ['hmmpanther','hmmpanther_start','hmmpanther_end'],
    ['blastprodom','blastprodom_start','blastprodom_end'],
    ['interpro','interpro_start','interpro_end'],
    ['signal_domain','signal_domain_start','signal_domain_end'],
    ['ncoils','ncoils_start','ncoils_end']
]

GENE_SCHEMA_COLUMNS = [
    'ensembl_gene_id',
    'entrez_gene',
    'gene_symbol',
    'chr',
    'strand',
    'genomic_start',
    'genomic_end'
]
GENE_SCHEMA = "(" + \
    "ensembl_gene_id text," + \
    "entrez_gene text," + \
    "gene_symbol text," + \
    "chr text," + \
    "strand integer," + \
    "genomic_start integer," + \
    "genomic_end integer" + \
    ")"

TRANSCRIPT_ANNOTATION_COLUMNS = [
    'ensembl_transcript_id',
    'ensembl_exon_id',
    'exon_chrom_start',
    'exon_chrom_end',
    'genomic_coding_start',
    'genomic_coding_end'
]
TRANSCRIPT_ANNOTATION_SCHEMA = "(" + \
    "ensembl_transcript_id text," + \
    "ensembl_exon_id text," + \
    "exon_chrom_start integer," + \
    "exon_chrom_end integer," + \
    "genomic_coding_start integer," + \
    "genomic_coding_end integer" + \
    ")"

TRANSCRIPT_SCHEMA_COLUMNS = [
    'ensembl_gene_id',
    'ensembl_transcript_id',
    'ensembl_protein_id',
    'transcript_biotype',
    'transcript_genomic_start',
    'transcript_genomic_end'
]
TRANSCRIPT_SCHEMA = "(" + \
    "ensembl_gene_id text," + \
    "ensembl_transcript_id text," + \
    "ensembl_protein_id text," + \
    "transcript_biotype text," + \
    "transcript_genomic_start integer," + \
    "transcript_genomic_end integer" + \
    ")"
