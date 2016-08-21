import multiprocessing
import sys
import os
import sqlite3
import requests

from tqdm import tqdm
from biomart import BiomartServer
from agfusion import utils

def split_into_n_lists(seq, n):
  avg = len(seq) / float(n)
  out = []
  last = 0.0

  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg

  return out

def _chunker(seq, n):
    """
    From a list return chunks of size n
    """

    return (seq[pos:pos + n] for pos in xrange(0, len(seq), n))

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
            line = line.decode('utf-8').split('\t')
            data.append(line)
            if len(line)!=7:
                print line
                print 'too short'


        return_dict[index]=data

    except requests.HTTPError:
        if ntries==3:
            print 'Max number of retries on this chunk!'
            raise
        else:
            print 'Chunk too large, splitting into two chunk and retrying...'
            for f in split_into_n_lists(filters[filters.keys()[0]],2):
                filter_sub = {filters.keys()[0]:f}
                _query_job(
                    index,
                    filter_sub,
                    attributes,
                    ntries+1,
                    ensembl,
                    return_dict
                )
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise

    return

class AGFusionSQlite3DB:
    """
    Class to handle methods around interacting with the AGFusion SQLite3
    database
    """

    def __init__(self,name,reference):

        self.name=name
        self.reference=reference
        self.datasets = {
            'GRCh38':'hsapiens_gene_ensembl',
            'GRCh37':'hsapiens_gene_ensembl',
            'GRCm38':'mmusculus_gene_ensembl'
        }

        if not os.path.exists(self.name):
            self.conn = sqlite3.connect(self.name)
            self.c = self.conn.cursor()
            self.initiate_tables()
        else:
            self.conn = sqlite3.connect(self.name)
            self.c = self.conn.cursor()

        self._ensembl=None
        self._biomart=None

    def _fetch(self,ids,filters,attributes,p,batch_size):
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

    def _fetch_gene_level_info(self,genes,p,table):

        print 'Fetching gene-level information...'

        data = self._fetch(
            ids=genes,
            filters='ensembl_gene_id',
            attributes=[
                'ensembl_gene_id',
                'entrezgene',
                'hgnc_symbol',
                'chromosome_name',
                'strand',
                'start_position',
                'end_position'
            ],
            p=p,
            batch_size=100
        )

        #process data

        print 'Adding gene-level information to the database...'

        #format the data correctly so it can be put into the database

        data_into_db = []
        for index, chunk in data.items():
            for r in chunk:
                data_into_db.append([
                    str(r[0]),
                    str(r[1]),
                    str(r[2]),
                    str(r[3]),
                    int(r[4]),
                    int(r[5]),
                    int(r[6]),
                ])

        self.c.execute('DELETE FROM ' + table)
        self.conn.commit()

        self.c.executemany('INSERT INTO ' + table + ' VALUES (?,?,?,?,?,?,?)', data_into_db)
        self.conn.commit()

    def _fetch_transcript_level_info(self,genes,p,table):

        print 'Fetching transcript-level information...'

        data = self._fetch(
            ids=genes,
            filters='ensembl_gene_id',
            attributes=[
                'ensembl_gene_id',
                'ensembl_transcript_id',
                'ensembl_peptide_id',
                'transcript_biotype',
                'transcript_start',
                'transcript_end',
                'transcript_length'
            ],
            p=p,
            batch_size=100
        )

        #process data

        print 'Adding transcript-level information to the database...'

        #format the data correctly so it can be put into the database

        data_into_db = []
        for index, chunk in data.items():
            for r in chunk:
                data_into_db.append([
                    str(r[0]),
                    str(r[1]),
                    str(r[2]),
                    str(r[3]),
                    int(r[4]),
                    int(r[5]),
                    int(r[6]),
                ])

        self.c.execute('DELETE FROM ' + table)
        self.conn.commit()

        self.c.executemany('INSERT INTO ' + table + ' VALUES (?,?,?,?,?,?,?)', data_into_db)
        self.conn.commit()

    def _fetch_protein_level_info(self,genes,p,table):
        print 'Fetching protein-level information...'

        data = self._fetch(
            ids=genes,
            filters='ensembl_transcript_id',
            attributes=[
                'ensembl_transcript_id',
                'ensembl_peptide_id',
                'transcript_biotype',
                'transcript_start',
                'transcript_end',
                'transcript_length'
            ],
            p=p,
            batch_size=100
        )

        #process data

        print 'Adding protein-level information to the database...'

        #format the data correctly so it can be put into the database

        data_into_db = []
        for index, chunk in data.items():
            for r in chunk:
                data_into_db.append([
                    str(r[0]),
                    str(r[1]),
                    str(r[2]),
                    str(r[3]),
                    int(r[4]),
                    int(r[5]),
                    int(r[6]),
                ])

        self.c.execute('DELETE FROM ' + table)
        self.conn.commit()

        self.c.executemany('INSERT INTO ' + table + ' VALUES (?,?,?,?,?,?,?)', data_into_db)
        self.conn.commit()

    def close(self):
        """
        Close the connection to the database
        """

        self.conn.close()

    def fetch_data(self,ensembl_server,p):

        self._biomart = BiomartServer(ensembl_server)
        self._ensembl = self._biomart.datasets[self.datasets[self.reference]]

        #get all the ensembl genes

        print 'Fetching all genes...'

        r=self._ensembl.search({'attributes':['ensembl_gene_id']})

        genes=[]

        for line in r.iter_lines():
            genes.append(line.decode('utf-8'))

        #self._fetch_gene_level_info(genes[0:1000],p,self.reference)

        self._fetch_transcript_level_info(genes[0:1000],p,self.reference + '_transcript')

        #get all transcript ids

        self.c.execute(
            "SELECT ensembl_transcript_id FROM " + \
            self.reference + '_transcript'
        )

        #self._fetch_protein_level_info(
        #    [str(i[0]) for i in self.c.fetchall()][0:1000],
        #    p,
        #    self.reference
        #)

    def initiate_tables(self):

        self.c.execute(
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
        self.c.execute(
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

        self.c.execute(
            "CREATE TABLE GRCh38_transcript (" + \
            "ensembl_gene_id text," + \
            "ensembl_transcript_id text," + \
            "ensembl_protein_id text," + \
            "transcript_biotype text," + \
            "transcript_genomic_start integer," + \
            "transcript_genomic_end integer," + \
            "transcript_length integer" + \
            ")"
        )
        self.c.execute(
            "CREATE TABLE GRCm38_transcript (" + \
            "ensembl_gene_id text," + \
            "ensembl_transcript_id text," + \
            "ensembl_protein_id text," + \
            "transcript_biotype text," + \
            "transcript_genomic_start integer," + \
            "transcript_genomic_end integer," + \
            "transcript_length integer" + \
            ")"
        )

        for i in utils.PROTEIN_DOMAIN:
            self.c.execute(
                "CREATE TABLE " + i[0] +" (" + \
                "ensembl_transcript_id text," +
                i[0] + " text," + \
                i[1] + " text," + \
                i[2] + " text" + \
                ")"
            )

        self.conn.commit()
