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

    def close(self):
        """
        Close the connection to the database
        """

        self.conn.close()
