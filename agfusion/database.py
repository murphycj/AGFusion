import multiprocessing
import sys
import os
import sqlite3

from biomart import BiomartServer
from agfusion import utils

def _chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

def _collect_results(result):
    results.extend(result)

def _query_job(filters,attributes,ntries):

    try:
        r=ensembl.search({
            'filters':filters,
            'attribute':attributes
        })

        data = []

        for line in r.iter_lines():
            data.append(line.decode('utf-8').split('\t'))
    except:
        if ntries==3:
            return None
        else:
            _query_job(filters,attributes,ntries+1)

    return data

class AGFusionSQlite3DB:
    """
    Class to handle methods around interacting with the AGFusion SQLite3
    database
    """

    def __init__(self,database_name):

        if not os.path.exists(database_name):
            print 'Database ' + database_name + ' does not exist!'
            sys.exit()

        self.database_name=database_name
        self.datasets = {
            'GRCh38':'hsapiens_gene_ensembl',
            'GRCh37':'hsapiens_gene_ensembl',
            'GRCm38':'mmusculus_gene_ensembl'
        }

        self.conn = sqlite3.connect(self.database_name)
        self.c = self.conn.cursor()

    def _fetch_gene_level_info(self,genes,p):
        attributes = [
            'ensembl_gene_id',
            'entrezgene',
            'hgnc_symbol'
            'chromosome_name',
            'chromosome_start',
            'chromosome_end',
            'strand'
        ]

        batch_size=100

        pool = multiprocessing.Pool(p)

        for i in _chunker(genes,batch_size):
            pool.apply_async(_query_job,args=({},[attributes],0,))

        pool.close()
        pool.join()

    def _fetch_transcript_level_info(self):
        attributes [
            'ensembl_gene_id',
            'ensembl_transcript_id',
            'ensembl_peptide_id',
        ]

    def _fetch_protein_level_info(self):
        pass

    def close(self):
        """
        Close the connection to the database
        """

        self.conn.close()

    def fetch_data(self,ensembl_server,ensembl_dataset,p):

        s = BiomartServer(ensembl_server)
        ensembl = s.datasets[ensembl_dataset]

        #get all the ensembl genes

        r=ensembl.search({'attributes':['ensembl_gene_id']})

        genes=[]

        for line in r.iter_lines():
            genes.append(line.decode('utf-8'))

        print len(genes)
        print 'Fetching gene-level information...'

        self._fetch_gene_level_info(genes,p)



    def initiate_tables(self):

        c.execute(
            "CREATE TABLE GRCh38 (" + \
            "ensembl_gene_id text," + \
            "entrez_gene text," + \
            "gene_symbol text," + \
            "chr text," + \
            "strand integer," + \
            "genomic_start integer," + \
            "genomic_end integer" + \
            ")"
        )
        c.execute(
            "CREATE TABLE GRCm38 (" + \
            "ensembl_gene_id text," + \
            "entrez_gene text," + \
            "gene_symbol text," + \
            "chr text," + \
            "strand integer," + \
            "genomic_start integer," + \
            "genomic_end integer" + \
            ")"
        )

        c.execute(
            "CREATE TABLE GRCh38_gene (" + \
            "ensembl_gene_id text," + \
            "ensembl_transcript_id text," + \
            "ensembl_protein_id text," + \
            "biotype text," + \
            "transcript_genomic_start integer," + \
            "transcript_genomic_end integer," + \
            "length integer" + \
            ")"
        )
        c.execute(
            "CREATE TABLE GRCm38_gene (" + \
            "ensembl_gene_id text," + \
            "ensembl_transcript_id text," + \
            "ensembl_protein_id text," + \
            "biotype text," + \
            "transcript_genomic_start integer," + \
            "transcript_genomic_end integer," + \
            "length integer" + \
            ")"
        )

        for i in agfusion.utils.PROTEIN_DOMAIN:
            c.execute(
                "CREATE TABLE " + i[0] +" (" + \
                "ensembl_transcript_id text," +
                i[0] + " text," + \
                i[1] + " text," + \
                i[2] + " text" + \
                ")"
            )

        self.conn.commit()
