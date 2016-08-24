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

PFAM_SCHEMA = "(" + \
    "ensembl_transcript_id text," + \
    "pfam text," + \
    "pfam_start text," + \
    "pfam_end text" + \
    ")"
PFAM_DOMAIN = ['pfam','pfam_start','pfam_end']

GENE_SCHEMA_COLUMNS = [
    'ensembl_gene_id',
    'entrez_gene',
    'gene_symbol'
]

GENE_SCHEMA = "(" + \
    "ensembl_gene_id text," + \
    "entrez_gene text," + \
    "gene_symbol text" + \
    ")"
