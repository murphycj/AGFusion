# all combinations except those that do not produce some sort of
# chimeric protein are commented out

import itertools

junction_locations = [
    'CDS','CDS (start)','CDS (end)','5UTR','5UTR (end)',
    '3UTR','3UTR (start)','exon'
    ]

CODING_COMBINATIONS = list(itertools.product(junction_locations,junction_locations))
CODING_COMBINATIONS = {i:False for i in CODING_COMBINATIONS}


CODING_COMBINATIONS['CDS','CDS']=True
CODING_COMBINATIONS['CDS','CDS (start)']=True
CODING_COMBINATIONS['CDS','CDS (end)']=True
CODING_COMBINATIONS['CDS','5UTR']=True
CODING_COMBINATIONS['CDS','5UTR (end)']=True
CODING_COMBINATIONS['CDS','3UTR']=True
CODING_COMBINATIONS['CDS','3UTR (start)']=True
CODING_COMBINATIONS['CDS (end)','CDS']=True
CODING_COMBINATIONS['CDS (end)','CDS (start)']=True
CODING_COMBINATIONS['CDS (end)','CDS (end)']=True
CODING_COMBINATIONS['CDS (end)','5UTR']=True
CODING_COMBINATIONS['CDS (end)','5UTR (end)']=True
CODING_COMBINATIONS['CDS (end)','3UTR']=True
CODING_COMBINATIONS['CDS (end)','3UTR (start)']=True
CODING_COMBINATIONS['CDS (start)','CDS (start)']=True
CODING_COMBINATIONS['5UTR','CDS (start)']=True
CODING_COMBINATIONS['5UTR (end)','CDS (start)']=True
CODING_COMBINATIONS['3UTR','CDS']=True
CODING_COMBINATIONS['3UTR','CDS (start)']=True
CODING_COMBINATIONS['3UTR','CDS (end)']=True
CODING_COMBINATIONS['3UTR','5UTR']=True
CODING_COMBINATIONS['3UTR','5UTR (end)']=True
CODING_COMBINATIONS['3UTR','3UTR']=True
CODING_COMBINATIONS['3UTR','3UTR (start)']=True
CODING_COMBINATIONS['3UTR (start)','CDS']=True
CODING_COMBINATIONS['3UTR (start)','CDS (start)']=True
CODING_COMBINATIONS['3UTR (start)','CDS (end)']=True
CODING_COMBINATIONS['3UTR (start)','5UTR']=True
CODING_COMBINATIONS['3UTR (start)','5UTR (end)']=True
CODING_COMBINATIONS['3UTR (start)','3UTR']=True
CODING_COMBINATIONS['3UTR (start)','3UTR (start)']=True


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
