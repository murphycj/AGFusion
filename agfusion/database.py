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

class AGFusionDB(object):
    """
    Class to handle methods around interacting with the AGFusion SQLite3
    database

    name: name of the database
    reference_name
    """

    def __init__(self,database_dir):

        self.database_dir=database_dir
        self.fastas = {}

        assert os.path.exists(os.path.abspath(os.path.join(self.database_dir,'agfusion.db'))), "database does not exist!"

        self.conn = sqlite3.connect(
            os.path.abspath(
                os.path.join(self.database_dir,'agfusion.db')
                )
            )
        self.c = self.conn.cursor()
        self.conn.commit()

    def _check_table(self,table):
        """
        Check if table exists
        """

        self.c.execute("SELECT * FROM sqlite_master WHERE name = \'" + table + "\' and type='table';")

        if len(self.c.fetchall())==0:
            return False
        else:
            return True

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

    def __init__(self,database_dir):

        if not os.path.exists(os.path.join(database_dir,'agfusion.db')):
            print 'Creating database...'
            fout = open(os.path.abspath(os.path.join(database_dir,'agfusion.db')),'a')
            fout.close()

        super(AGFusionDBBManager,self).__init__(database_dir)

        self._ensembl=None
        self._biomart=None
        #self.reference=reference

        self._check_for_tables()

    def _check_for_tables(self):

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

    def _fetch_protein_level_info(self,genes,p,columns,table):
        print 'Fetching protein-level information...'

        data = self._fetch(
            ids=genes,
            filters='ensembl_transcript_id',
            attributes=['ensembl_transcript_id'] + columns,
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
                    str(r[3])
                ])

        self.c.execute('DELETE FROM ' + table)
        self.conn.commit()

        self.c.executemany('INSERT INTO ' + table + ' VALUES (?,?,?,?)', data_into_db)
        self.conn.commit()

    def fetch_data(self,ensembl_server,ensembl_dataset,p,transcripts):

        self._biomart = BiomartServer(ensembl_server)
        self._ensembl = self._biomart.datasets[ensembl_dataset]

        self._fetch_protein_level_info(
            transcripts,
            p,
            utils.PFAM_DOMAIN,
            'pfam'
        )
