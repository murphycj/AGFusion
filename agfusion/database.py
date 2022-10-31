"""Classes for AGFusion database objects
"""
import logging
import sqlite3
import sys
from os.path import abspath, exists, join

from agfusion import utils

logger = logging.getLogger("AGFusion")

try:
    import MySQLdb
except ImportError:
    logger.debug("Could not import MySQLdb.")


class AGFusionDB:
    """
    Class to handle methods around interacting with the AGFusion SQLite3
    database

    name: name of the database
    reference_name
    """

    def __init__(self, database=None, debug=False):

        self.database = abspath(database)
        self.fastas = {}

        self.logger = logging.getLogger("AGFusion")
        if debug:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        if debug:
            ch.setLevel(logging.DEBUG)
        else:
            ch.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        assert exists(self.database), (
            f"AGFusion database at {database} does not exist! "
            "Either run 'agfusion download' or specify the location of the AGFusion "
            "database with the --dbpath flag."
        )

        self.sqlite3_db = sqlite3.connect(abspath(self.database))
        self.sqlite3_cursor = self.sqlite3_db.cursor()
        self.sqlite3_db.commit()

        self.logger.debug("Connected to the database %s", abspath(self.database))

        self.build = ""


class AGFusionDBBManager:
    """
    Class to handle methods managing the SQLite3 database

    name: name of the database
    reference
    """

    def __init__(self, db_dir, species, release, pfam, server):

        self.species = species
        self.release = release
        self.server = server
        self.build = self.species + "_" + str(self.release)
        self.table = utils.ENSEMBL_MYSQL_TABLES[self.species][self.release]

        self.database = join(
            abspath(db_dir),
            "agfusion." + self.species + "." + str(self.release) + ".db",
        )
        self.fastas = {}

        self.logger = logging.getLogger("AGFusion")
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        if not exists(self.database):
            fout = open(abspath(self.database), "a")
            fout.close()
            self.logger.info("Created the database %s", abspath(self.database))

        self.sqlite3_db = sqlite3.connect(abspath(self.database))
        self.sqlite3_cursor = self.sqlite3_db.cursor()
        self.sqlite3_db.commit()

        self.logger.info("Connected to the database %s", abspath(self.database))

        self.ensembl_db = MySQLdb.connect(server, "anonymous")
        self.ensembl_cursor = self.ensembl_db.cursor()
        self.ensembl_cursor.execute("use " + self.table + ";")

        self.logger.info("Connected to the ensembl MySQL server at %s", server)
        self.logger.info("MySQL - use %s", self.table + ";")

        self._check_for_tables()

        # load the PFAM mappings

        self.pfam_mapping = {}

        for line in open(pfam, "r"):
            line = line.rstrip().split("\t")

            pfam_id = line[0]
            pfam_name = line[1]
            pfam_desc = line[3]

            self.pfam_mapping[pfam_id] = {"name": pfam_name, "desc": pfam_desc}

    def _check_for_tables(self):
        """Check if table exists."""

        self.sqlite3_cursor.execute("drop table if exists " + self.build)

        sqlite3_command = (
            "CREATE TABLE "
            + self.build
            + " ("
            + "gene_id text,"
            + "stable_id text,"
            + "entrez_id text,"
            + "gene_name text,"
            + "canonical_transcript_id text);"
        )

        self.logger.info("MySQL - %s", sqlite3_command)

        self.sqlite3_cursor.execute(sqlite3_command)
        self.sqlite3_db.commit()

        # transcript table

        self.sqlite3_cursor.execute("drop table if exists " + self.build + "_transcript")

        sqlite3_command = (
            "CREATE TABLE "
            + self.build
            + "_transcript ("
            + "transcript_id text,"
            + "gene_id text,"
            + "transcript_stable_id text,"
            + "translation_id text);"
        )

        self.logger.info("MySQL - %s", sqlite3_command)

        self.sqlite3_cursor.execute(sqlite3_command)
        self.sqlite3_db.commit()

        # refseq table

        self.sqlite3_cursor.execute("drop table if exists " + self.build + "_refseq")

        sqlite3_command = (
            "CREATE TABLE "
            + self.build
            + "_refseq ("
            + "transcript_id text,"
            + "transcript_stable_id text,"
            + "refseq_id text);"
        )

        self.logger.info("MySQL - %s", sqlite3_command)

        self.sqlite3_cursor.execute(sqlite3_command)
        self.sqlite3_db.commit()

        # protein annotation tables

        for protein_annotation in utils.PROTEIN_ANNOTATIONS:
            self.sqlite3_cursor.execute(
                "drop table if exists " + self.build + "_" + protein_annotation
            )

            sqlite3_command = (
                f"CREATE TABLE {self.build}_{protein_annotation} ("
                "translation_id text,stable_id text,hit_id text,"
                "seq_start integer,seq_end integer,hit_description text,"
                "hit_name text);"
            )

            self.logger.info("MySQL - %s", sqlite3_command)

            self.sqlite3_cursor.execute(sqlite3_command)
            self.sqlite3_db.commit()

    def fetch_gene_names(self):
        """
        Fetch gene by name.
        """

        genes = {}

        # fetch all gene stable ids

        if self.release < 65:
            mysql_command = (
                "SELECT gene.gene_id, gene_stable_id.stable_id, gene.canonical_transcript_id "
                "FROM gene, gene_stable_id WHERE gene.gene_id = gene_stable_id.gene_id;"
            )
        else:
            mysql_command = (
                "SELECT gene.gene_id, gene.stable_id, gene.canonical_transcript_id FROM gene;"
            )

        self.logger.info("MySQL - %s", mysql_command)

        self.ensembl_cursor.execute(mysql_command)
        for g in self.ensembl_cursor.fetchall():
            genes[g[0]] = {
                "stable_id": g[1],
                "canonical_transcript_id": g[2],
                "entrez_id": "",
                "gene_name": "",
            }

        # fetch entrez IDS

        mysql_command = (
            "SELECT gene.gene_id, xref.dbprimary_acc FROM "
            "gene, object_xref, xref,external_db WHERE "
            "gene.gene_id = object_xref.ensembl_id AND "
            "object_xref.ensembl_object_type = 'Gene' AND "
            "object_xref.xref_id = xref.xref_id AND "
            "xref.external_db_id = external_db.external_db_id AND "
            "external_db.db_name = 'EntrezGene';"
        )
        self.logger.info("MySQL - %s", mysql_command)

        self.ensembl_cursor.execute(mysql_command)

        for g in self.ensembl_cursor.fetchall():
            genes[g[0]]["entrez_id"] = g[1]

        # fetch gene names

        if self.build.find("homo_sapiens") != -1:
            gene_name_db = "HGNC"
        elif self.build.find("mus_musculus") != -1:
            gene_name_db = "MGI"
        else:
            self.logger.error("Do not know where to fetch gene names for %s", self.build)
            sys.exit(1)

        mysql_command = (
            """SELECT gene.gene_id, xref.display_label FROM " \
            "gene, object_xref, xref,external_db WHERE " \
            "gene.gene_id = object_xref.ensembl_id AND " \
            "object_xref.ensembl_object_type = 'Gene' AND " \
            "object_xref.xref_id = xref.xref_id AND " \
            "xref.external_db_id = external_db.external_db_id AND " \
            "external_db.db_name = '"""
            + gene_name_db
            + """';"""
        )

        self.logger.info("MySQL - %s", mysql_command)

        self.ensembl_cursor.execute(mysql_command)

        for g, e in self.ensembl_cursor.fetchall():
            genes[g]["gene_name"] = e

        genes = [
            [
                i,
                genes[i]["stable_id"],
                genes[i]["entrez_id"],
                genes[i]["gene_name"],
                genes[i]["canonical_transcript_id"],
            ]
            for i in list(genes.keys())
        ]

        command = (
            f"SQLite - INSERT INTO {self.build} VALUES "
            "(gene_id,stable_id,entrez_id,gene_name,canonical_transcript_id)"
        )
        self.logger.info(command)

        self.sqlite3_cursor.executemany("INSERT INTO " + self.build + " VALUES (?,?,?,?,?)", genes)

        self.sqlite3_db.commit()

    def fetch_transcript_table(self):
        """
        Fetch all transcripts
        """

        if self.release < 65:
            mysql_command = (
                "SELECT transcript.transcript_id, transcript.gene_id, "
                "transcript_stable_id.stable_id FROM transcript, "
                "transcript_stable_id WHERE "
                "transcript.transcript_id = transcript_stable_id.transcript_id;"
            )
        else:
            mysql_command = (
                "SELECT transcript.transcript_id, "
                "transcript.gene_id, transcript.stable_id FROM transcript;"
            )

        self.logger.info("MySQL - %s", mysql_command)
        self.ensembl_cursor.execute(mysql_command)
        transcripts = {}
        for transcript in self.ensembl_cursor.fetchall():
            transcripts[transcript[0]] = {
                "gene_id": transcript[1],
                "transcript_stable_id": transcript[2],
                "translation_id": "",
            }

        # Fetch all transcripts with tranlstions

        mysql_command = (
            "SELECT transcript.transcript_id, "
            "translation.translation_id FROM "
            "transcript, translation WHERE "
            "transcript.transcript_id = translation.transcript_id;"
        )
        self.logger.info("MySQL - %s", mysql_command)
        self.ensembl_cursor.execute(mysql_command)
        for transcript in self.ensembl_cursor.fetchall():
            transcripts[transcript[0]]["translation_id"] = transcript[1]

        transcripts = [
            [
                i,
                transcripts[i]["gene_id"],
                transcripts[i]["transcript_stable_id"],
                transcripts[i]["translation_id"],
            ]
            for i in list(transcripts.keys())
        ]

        log = (
            f"SQLite - INSERT INTO {self.build}_transcript "
            "VALUES (transcript_id,gene_id,transcript_stable_id,translation_id)"
        )
        self.logger.info(log)

        self.sqlite3_cursor.executemany(
            "INSERT INTO " + self.build + "_transcript VALUES (?,?,?,?)", transcripts
        )

        self.sqlite3_db.commit()

    def fetch_refseq_table(self):
        """
        fetch RefSeq IDS
        """

        if self.release < 65:
            mysql_command = (
                "SELECT transcript.transcript_id, transcript_stable_id.stable_id, "
                "xref.display_label FROM "
                "transcript, transcript_stable_id, object_xref, xref, external_db WHERE "
                "transcript.transcript_id = transcript_stable_id.transcript_id and "
                "transcript.transcript_id = object_xref.ensembl_id AND "
                "object_xref.ensembl_object_type = 'Transcript' AND "
                "object_xref.xref_id = xref.xref_id AND "
                "xref.external_db_id = external_db.external_db_id AND "
                "external_db.db_name = 'RefSeq_mRNA';"
            )
        else:
            mysql_command = (
                "SELECT transcript.transcript_id, transcript.stable_id, xref.display_label "
                "FROM transcript, object_xref, xref, external_db WHERE "
                "transcript.transcript_id = object_xref.ensembl_id AND "
                "object_xref.ensembl_object_type = 'Transcript' AND "
                "object_xref.xref_id = xref.xref_id AND "
                "xref.external_db_id = external_db.external_db_id AND "
                "external_db.db_name = 'RefSeq_mRNA';"
            )
        self.logger.info("MySQL - %s", mysql_command)
        self.ensembl_cursor.execute(mysql_command)

        refseqs = [[i[0], i[1], i[2]] for i in self.ensembl_cursor.fetchall()]

        log = (
            f"SQLite - INSERT INTO {self.build}_refseq VALUES "
            "(transcript_id,transcript_stable_id,refseq_id)"
        )
        self.logger.info(log)

        self.sqlite3_cursor.executemany(
            "INSERT INTO " + self.build + "_refseq VALUES (?,?,?)", refseqs
        )

        self.sqlite3_db.commit()

    def fetch_protein_annotation(self):
        """Fetch protein annotation information."""

        for protein_annotation in utils.PROTEIN_ANNOTATIONS:

            if self.release < 70:
                mysql_command = (
                    "SELECT translation.translation_id, translation_stable_id.stable_id, "
                    "protein_feature.hit_name, protein_feature.seq_start, "
                    "protein_feature.seq_end FROM "
                    "analysis, analysis_description, protein_feature, "
                    "translation, translation_stable_id WHERE "
                    "translation_stable_id.translation_id = translation.translation_id AND "
                    "protein_feature.translation_id = translation.translation_id AND "
                    "protein_feature.analysis_id = analysis.analysis_id AND "
                    "analysis.analysis_id = analysis_description.analysis_id AND "
                    "analysis.logic_name = '" + protein_annotation + "';"
                )
            else:
                mysql_command = (
                    "SELECT translation.translation_id, translation.stable_id, "
                    "protein_feature.hit_name, protein_feature.seq_start, "
                    "protein_feature.seq_end, protein_feature.hit_description FROM "
                    "analysis, analysis_description, protein_feature, translation WHERE "
                    "protein_feature.translation_id = translation.translation_id AND "
                    "protein_feature.analysis_id = analysis.analysis_id AND "
                    "analysis.analysis_id = analysis_description.analysis_id AND "
                    "analysis.logic_name = '" + protein_annotation + "';"
                )

            self.logger.info("MySQL - %s", mysql_command)

            self.ensembl_cursor.execute(mysql_command)

            data = list(self.ensembl_cursor.fetchall())

            if self.release < 70:
                data = [i + (None,) for i in data]

            if protein_annotation == "pfam":
                for i, tmp in enumerate(data):
                    tmp = list(tmp)
                    if tmp[2] in self.pfam_mapping:
                        tmp.append(self.pfam_mapping[tmp[2]]["name"])
                        tmp[5] = self.pfam_mapping[tmp[2]]["desc"]
                    else:
                        tmp.append("")
                    data[i] = tmp
            else:
                for i, tmp in enumerate(data):
                    tmp = list(tmp)
                    tmp.append(None)
                    data[i] = tmp

            msg = (
                f"SQLite - INSERT INTO {self.build}_{protein_annotation}"
                " VALUES "
                "(translation_id,stable_id,hit_id,seq_start,seq_end,hit_description,hit_name)"
            )

            self.logger.info(msg)

            self.sqlite3_cursor.executemany(
                "INSERT INTO " + self.build + "_" + protein_annotation + " VALUES (?,?,?,?,?,?,?)",
                data,
            )

            self.sqlite3_db.commit()
