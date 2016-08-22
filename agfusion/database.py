import multiprocessing
import sys
import os
import sqlite3
import requests
import glob

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
            #if len(line)!=7:
                #print line
                #print 'too short'


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

class AGFusionDB:
    """
    Class to handle methods around interacting with the AGFusion SQLite3
    database

    name: name of the database
    reference_name
    """

    def __init__(self,database_name):

        self.database_name=database_name
        self.fastas = {}

        self.conn = sqlite3.connect(self.database_name)
        self.c = self.conn.cursor()
        self.conn.commit()

        #check that the files that are referenced in the DB are where
        #they are supposed to be

        self.c.execute('SELECT * FROM FastaGtf')
        files = self.c.fetchall()
        for ref in files:
            for f in ref:
                assert os.path.exists(f), "The file %s cannot be found!" % f


    def close(self):
        """
        Close the connection to the database
        """

        self.conn.close()

class AGFusionDBBManager(AGFusionDB):
    """
    Class to handle methods managing the SQLite3 database

    name: name of the database
    reference
    """

    def __init__(self,database_name):

        if not os.path.exists(self.database_name):
            print 'Creating database...'

        super(AGFusionDB,self).__init__(database_name)

        self._ensembl=None
        self._biomart=None

    def _check_table(self,table):
        """
        Check if table exists
        """

        self.c.execute("SELECT * FROM sqlite_master WHERE name = \'" + table + "\' and type='table';")

        if len(db.c.fetchall())==0:
            return False
        else:
            return True

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

    def _fetch_transcript_level_info(self,genes,p):

        print 'Fetching transcript-level information...'

        data = self._fetch(
            ids=genes,
            filters='ensembl_gene_id',
            attributes=[
                'ensembl_gene_id',
                'ensembl_transcript_id',
                'ensembl_peptide_id',
                'ensembl_exon_id',
                'transcript_biotype',
                'transcript_start',
                'transcript_end',
                'exon_chrom_start',
                'exon_chrom_end',
                'genomic_coding_start',
                'genomic_coding_end'
            ],
            p=p,
            batch_size=100
        )

        #process data

        print 'Adding transcript-level information to the database...'

        #insert into transcipt table

        #format the data correctly so it can be put into the database

        data_into_db = []
        for index, chunk in data.items():
            for r in chunk:
                temp = [
                    str(r[0]),
                    str(r[1]),
                    str(r[2]),
                    str(r[4]),
                    -1 if r[5]==u'' else int(r[5]),
                    -1 if r[6]==u'' else int(r[6]),
                ]
                if temp not in data_into_db:
                    data_into_db.append(temp)

        self.c.execute('DELETE FROM ' + self.reference + '_transcript')
        self.conn.commit()

        self.c.executemany(
            'INSERT INTO ' + self.reference + '_transcript VALUES (?,?,?,?,?,?)',
            data_into_db
        )
        self.conn.commit()

        #insert into transcipt_annotation table

        #format the data correctly so it can be put into the database

        data_into_db = []
        for index, chunk in data.items():
            for r in chunk:
                data_into_db.append([
                    str(r[1]),
                    str(r[3]),
                    -1 if r[7]==u'' else int(r[7]),
                    -1 if r[8]==u'' else int(r[8]),
                    -1 if r[9]==u'' else int(r[9]),
                    -1 if r[10]==u'' else int(r[10]),
                ])

        self.c.execute('DELETE FROM ' + self.reference + '_annotation_transcript')
        self.conn.commit()

        self.c.executemany(
            'INSERT INTO ' + self.reference + '_annotation_transcript VALUES (?,?,?,?,?,?)',
            data_into_db
        )
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

    def fetch_data(self,ensembl_server,ensembl_dataset,reference,reference_dir,p,reference):

        self._check_for_tables(reference,reference_dir)

        self._biomart = BiomartServer(ensembl_server)
        self._ensembl = self._biomart.datasets[ensembl_dataset]

        #get all the ensembl genes

        #print 'Fetching all genes...'

        #r=self._ensembl.search({'attributes':['ensembl_gene_id']})

        #genes=[]

        #for line in r.iter_lines():
        #    genes.append(line.decode('utf-8'))

        #self._fetch_gene_level_info(["ENSG00000066468","ENSG00000197959"],p,reference)
        #self._fetch_gene_level_info(genes[0:1000],p,self.reference)

        #self._fetch_transcript_level_info(["ENSG00000066468","ENSG00000197959"],p)
        #self._fetch_transcript_level_info(genes,p)

        #get all transcript ids

        self.c.execute(
            "SELECT ensembl_transcript_id FROM " + \
            reference + '_transcript'
        )

        #self._fetch_protein_level_info(
        #    [str(i[0]) for i in self.c.fetchall()][0:1000],
        #    p,
        #    self.reference
        #)

    def add_fasta_gtf(self,reference,reference_dir):
        """
        Locate the fasta and gtf files and add the path to them to the database
        """

        file_types = ['primary_assembly.fa','cdna.all.fa','cds.all.fa','pep.all.fa','gtf']
        files = []
        for f in file_types:
            r = glob.glob(os.path.join(reference_dir,reference + '*' * f))

            assert len(r)==1, "No files or too many reference files were found %s" % r

            files += r

        self.c.executemany(
            "INSERT INTO FastaGtf VALUES (?,?,?,?,?)",
            files
        )
        self.conn.commit()

    def chack_tables(self,reference,reference_dir):

        if not self._check_table(reference_name):
            self.c.execute(
                "CREATE TABLE " + reference_name + " " + utils.GENE_SCHEMA
            )
            self.conn.commit()

        if not self._check_table(reference_name + '_transcript'):
            self.c.execute(
                "CREATE TABLE " + reference_name + "_transcript " + utils.TRANSCRIPT_SCHEMA
            )
            self.conn.commit()

        if not self._check_table(reference_name + '_annotation_transcript'):
            self.c.execute(
                "CREATE TABLE " + reference_name + "_annotation_transcript " + utils.TRANSCRIPT_ANNOTATION_SCHEMA
            )
            self.conn.commit()


        for i in utils.PROTEIN_DOMAIN:
            if not self._check_table(i[0]):
                self.c.execute(
                    "CREATE TABLE " + i[0] +" (" + \
                    "ensembl_transcript_id text," +
                    i[0] + " text," + \
                    i[1] + " text," + \
                    i[2] + " text" + \
                    ")"
                )
                self.conn.commit()
