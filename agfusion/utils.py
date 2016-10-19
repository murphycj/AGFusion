# all combinations except those that do not produce some sort of
# chimeric protein are commented out

import itertools

junction_locations = [
    'CDS','CDS (start)','CDS (end)','5UTR','5UTR (end)',
    '3UTR','3UTR (start)','exon','intron','intron (cds)',
    'intron (before cds)','intron (after cds)'
    ]

CODING_COMBINATIONS = list(itertools.product(junction_locations,junction_locations))
CODING_COMBINATIONS = {i:{'protein_coding_potential':False,'truncated':False} for i in CODING_COMBINATIONS}


CODING_COMBINATIONS['CDS','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS','CDS (end)'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS','5UTR'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS','5UTR (end)'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS','3UTR'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS','3UTR (start)'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS','intron'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS','intron (before cds)'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS','intron (after cds)'] = {'protein_coding_potential':True,'truncated':True}

CODING_COMBINATIONS['CDS (end)','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','5UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','3UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','intron'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (end)','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}

#CODING_COMBINATIONS['CDS (start)','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (start)','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['CDS (start)','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['CDS (start)','5UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['CDS (start)','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['CDS (start)','3UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['CDS (start)','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (start)','intron'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS (start)','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['CDS (start)','intron (before cds)'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['CDS (start)','intron (after cds)'] = {'protein_coding_potential':True,'truncated':True}

#CODING_COMBINATIONS['5UTR','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['5UTR','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR','5UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR','3UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR','intron'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['5UTR','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}

#CODING_COMBINATIONS['5UTR (end)','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['5UTR (end)','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR (end)','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR (end)','5UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR (end)','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR (end)','3UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR (end)','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR (end)','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['5UTR (end)','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['5UTR (end)','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}

CODING_COMBINATIONS['3UTR','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','5UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','3UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','intron'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}

CODING_COMBINATIONS['3UTR (start)','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','5UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','3UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','intron'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (start)','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}

#CODING_COMBINATIONS['3UTR (end)','CDS'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','5UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','3UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','intron'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['3UTR (end)','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['3UTR (end)','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}

#CODING_COMBINATIONS['intron','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron','5UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron','3UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron','intron'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron','intron (after cds)'] = {'protein_coding_potential':{'protein_coding_potential':True,'truncated':False},'truncated':False}

CODING_COMBINATIONS['intron (cds)','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (cds)','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (cds)','CDS (end)'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['intron (cds)','5UTR'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['intron (cds)','5UTR (end)'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['intron (cds)','3UTR'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['intron (cds)','3UTR (start)'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['intron (cds)','intron'] = {'protein_coding_potential':True,'truncated':True}
CODING_COMBINATIONS['intron (cds)','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron (cds)','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron (cds)','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}

#CODING_COMBINATIONS['intron (before cds)','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (before cds)','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron (before cds)','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (before cds)','5UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (before cds)','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron (before cds)','3UTR'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron (before cds)','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron (before cds)','intron'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron (before cds)','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (before cds)','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
#CODING_COMBINATIONS['intron (before cds)','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}

CODING_COMBINATIONS['intron (after cds)','CDS'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','CDS (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','CDS (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','5UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','5UTR (end)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','3UTR'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','3UTR (start)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','intron'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','intron (cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','intron (before cds)'] = {'protein_coding_potential':True,'truncated':False}
CODING_COMBINATIONS['intron (after cds)','intron (after cds)'] = {'protein_coding_potential':True,'truncated':False}

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
