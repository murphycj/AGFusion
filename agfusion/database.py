import sys
import os
import sqlite3
import requests
import logging


class AGFusionDB(object):
    """
    Class to handle methods around interacting with the AGFusion SQLite3
    database

    name: name of the database
    reference_name
    """

    def __init__(self, database):

        self.database = os.path.abspath(database)
        self.fastas = {}

        self.logger = logging.getLogger('AGFusion')
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        assert os.path.exists(self.database), "database does not exist!"

        #self.conn = sqlite3.connect(
        #    self.database
        #)
        #self.c = self.conn.cursor()
        #self.conn.commit()

        self.sqlite3_db = sqlite3.connect(
            os.path.abspath(self.database)
        )
        self.sqlite3_cursor = self.sqlite3_db.cursor()
        self.sqlite3_db.commit()

        self.logger.info(
            'Connected to the database ' + os.path.abspath(self.database)
        )

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

    def __init__(self, database, build, server):

        super(AGFusionDBBManager, self).__init__(database)
        self.build = build
        self.server = server

        if not os.path.exists(self.database):
            fout = open(os.path.abspath(self.database), 'a')
            fout.close()
            self.logger.info(
                "Created the database " + os.path.abspath(self.database)
            )

        import MySQLdb

        self.ensembl_db = MySQLdb.connect(server, 'anonymous')
        self.ensembl_cursor = self.ensembl_db.cursor()
        self.ensembl_cursor.execute('use ' + build + ';')

        self.logger.info(
            'Connected to the ensembl MySQL server at ' + server
        )
        self.logger.info('MySQL - use ' + build + ';')

        self._check_for_tables()

    def _check_for_tables(self):

        self.sqlite3_cursor.execute('drop table if exists ' + self.build)

        sqlite3_command = "CREATE TABLE " + self.build + " (" + \
            "gene_id text," + \
            "stable_id text," + \
            "entrez_id text," + \
            "gene_name text," + \
            "canonical_transcript_id text);"

        self.logger.info('SQLite - ' + sqlite3_command)

        self.sqlite3_cursor.execute(
            sqlite3_command
        )
        self.sqlite3_db.commit()

        #transcript table

        self.sqlite3_cursor.execute('drop table if exists ' + self.build + '_transcript')

        sqlite3_command = "CREATE TABLE " + self.build + "_transcript (" + \
            "transcript_id text," + \
            "transcript_stable_id text," + \
            "refseq_id text," + \
            "translation_id text);"

        self.logger.info('SQLite - ' + sqlite3_command)

        self.sqlite3_cursor.execute(
            sqlite3_command
        )
        self.sqlite3_db.commit()

    def fetch_refseq(self):
        pass

    def fetch_gene_names(self):

        genes = {}

        # fetch all gene stable ids

        mysql_command = "SELECT gene.gene_id, gene.stable_id, " + \
            "gene.canonical_transcript_id FROM gene;"

        self.logger.info('MySQL - ' + mysql_command)

        self.ensembl_cursor.execute(
            mysql_command
        )
        for g in self.ensembl_cursor.fetchall():
            genes[g[0]] = {
                'stable_id': g[1],
                'canonical_transcript_id': g[2],
                'entrez_id': '',
                'gene_name': ''
            }

        # fetch entrez IDS

        mysql_command = "SELECT gene.gene_id, xref.dbprimary_acc " + \
            "FROM gene, object_xref, xref,external_db " + \
            "WHERE gene.gene_id = object_xref.ensembl_id " + \
            "AND object_xref.ensembl_object_type = 'Gene' " + \
            "AND object_xref.xref_id = xref.xref_id " + \
            "AND xref.external_db_id = external_db.external_db_id " + \
            "AND external_db.db_name = 'EntrezGene';"

        self.logger.info('MySQL - ' + mysql_command)

        self.ensembl_cursor.execute(mysql_command)

        for g in self.ensembl_cursor.fetchall():
            genes[g[0]]['entrez_id'] = g[1]

        # fetch gene names

        if self.build.find('homo_sapiens') != -1:
            gene_name_db = 'HGNC'
        elif self.build.find('mus_musculus') != -1:
            gene_name_db = 'MGI'
        else:
            self.logger.error(
                'Do not know where to fetch gene names for ' + self.build
            )
            sys.exit()

        mysql_command = \
            """ \
            SELECT gene.gene_id, xref.display_label \
            FROM gene, object_xref, xref,external_db \
            WHERE gene.gene_id = object_xref.ensembl_id \
            AND object_xref.ensembl_object_type = 'Gene' \
            AND object_xref.xref_id = xref.xref_id \
            AND xref.external_db_id = external_db.external_db_id \
            AND external_db.db_name = '""" + gene_name_db + """'; \
            """

        self.logger.info('MySQL - ' + mysql_command)

        self.ensembl_cursor.execute(mysql_command)

        for g, e in self.ensembl_cursor.fetchall():
            genes[g]['gene_name'] = e

        genes = [
            [
                i,
                genes[i]['stable_id'],
                genes[i]['entrez_id'],
                genes[i]['gene_name'],
                genes[i]['canonical_transcript_id']
            ] for i in genes.keys()
        ]

        self.logger.info(
            'SQLite - INSERT INTO ' +
            self.build + ' VALUES (gene_id,stable_id,entrez_id,gene_name,canonical_transcript_id)'
        )

        self.sqlite3_cursor.executemany(
            'INSERT INTO ' + self.build + ' VALUES (?,?,?,?,?)', genes
        )

        self.sqlite3_db.commit()

    def fetch_transcript_table(self):
        mysql_command = "SELECT transcript.transcript_id, transcript.stable_id, translation.translation_id FROM transcript, translation WHERE transcript.transcript_id = translation.transcript_id;"

        self.logger.info('MySQL - ' + mysql_command)

        self.ensembl_cursor.execute(
            mysql_command
        )
        transcripts = {}
        for transcript in self.ensembl_cursor.fetchall():
            transcripts[transcript[0]] = {
                'transcript_stable_id': transcript[1],
                'translation_id': transcript[2],
                'refseq_id': ''
            }

        mysql_command = "SELECT transcript.transcript_id, xref.display_label FROM transcript, object_xref, xref, external_db WHERE transcript.transcript_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = \'Transcript\' AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id AND external_db.db_name = \'RefSeq_mRNA\';"
        self.logger.info('MySQL - ' + mysql_command)
        self.ensembl_cursor.execute(
            mysql_command
        )

        for transcript in self.ensembl_cursor.fetchall():
            transcripts[transcript[0]]['refseq_id'] = transcript[1]

        transcripts = [
            [
                i,
                transcripts[i]['transcript_stable_id'],
                transcripts[i]['translation_id'],
                transcripts[i]['refseq_id']
            ] for i in transcripts.keys()
        ]

        self.logger.info(
            'SQLite - INSERT INTO ' +
            self.build + '_transcript VALUES (transcript_id,transcript_stable_id,refseq_id,translation_id)'
        )

        self.sqlite3_cursor.executemany(
            'INSERT INTO ' + self.build + '_transcript VALUES (?,?,?,?)', transcripts
        )

        self.sqlite3_db.commit()


    def fetch_protein_annotation(self):

        for protein_annotation in ['PFAM','SMART','Superfamily','TIGR','HAMAP','ProSite']
