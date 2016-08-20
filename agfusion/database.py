import multiprocessing
import sys
import os
import sqlite3

from tqdm import tqdm
from biomart import BiomartServer
from agfusion import utils

def _chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

def _collect_results(result):

    if result is None:
        print 'Could not fetch some data for some reason...'
        sys.exit()
        
    results.extend(result)

def _query_job(index,filters,attributes,ntries,ensembl,return_dict):

    data = []

    try:

        r=ensembl.search({
            'filters':filters,
            'attributes':attributes[0]
        })

        for line in r.iter_lines():
            data.append(line.decode('utf-8').split('\t'))

        return_dict[index]=data

    except:
        print 'Retrying...'
        if ntries==3:
            return None
        else:
            _query_job(index,filters,attributes,ntries+1,ensembl,return_dict)

    return

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

        self._ensembl=None
        self._biomart=None

    def _fetch(self,ids,filters,attrubutes,p,batch_size):
        """
        Abstract method for fetching data in batches from the biomart server.
        The data is fetched in batches because the biomart python package does
        not seem to like having to fetch data in too large of chunks
        """

        #setup multiprocessing

        pool = multiprocessing.Pool(p)
        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        #chunk the data

        chunks = [i for i in _chunker(ids,batch_size)]
        sub_chunks = [i for i in _chunker(chunks,p)]

        index=0

        #fetch the data

        for chunk in tqdm(sub_chunks):
            jobs = []
            for id_set in chunk:
                p = multiprocessing.Process(target=_query_job,args=(index,{filters:id_set},[attributes],0,self._ensembl,return_dict,))
                jobs.append(p)
                p.start()
                index+=1

            for j in jobs:
                j.join()

        return return_dict


    def _fetch_gene_level_info(self,genes,p):

        print 'Fetching gene-level information...'

        data = self._fetch(
            ids=genes,
            filters='ensembl_gene_id',
            attrubutes=[
                'ensembl_gene_id',
                'entrezgene',
                'hgnc_symbol',
                'chromosome_name',
                'start_position',
                'end_position',
                'strand'
            ],
            p=p,
            batch_size=100
        )

        #process data

        print 'Adding gene-level information to the database...'

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

        self._biomart = BiomartServer(ensembl_server)
        self._ensembl = self._biomart.datasets[ensembl_dataset]

        #get all the ensembl genes

        print 'Fetching all genes...'

        r=self._ensembl.search({'attributes':['ensembl_gene_id']})

        genes=[]

        for line in r.iter_lines():
            genes.append(line.decode('utf-8'))

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
