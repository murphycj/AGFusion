import sys
import os
import sqlite3
import requests

class AGFusionDB(object):
    """
    Class to handle methods around interacting with the AGFusion SQLite3
    database

    name: name of the database
    reference_name
    """

    def __init__(self,database):

        self.database=os.path.abspath(database)
        self.fastas = {}

        self.logger = logging.getLogger('AGFusion')
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        assert os.path.exists(self.database), "database does not exist!"

        self.conn = sqlite3.connect(
            self.database
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
