import sqlite3
import os
import sys
import argparse

from agfusion import database, model

def build_database():

    parser = argparse.ArgumentParser(
        description='Build or update the SQLite3 database for human and mouse ' + \
        'genomes by querying Biomart. If the database given by --name already ' + \
        'exists, then species-specific portion of the database will be updated.'
    )
    parser.add_argument(
        '--name',
        type=str,
        required=True,
        help='Name of the database (default: updates the current database if it already exists)'
    )
    parser.add_argument(
        '--p',
        type=int,
        default=1,
        help='Number of processers to use to fetch data from Biomart (default 1).'
    )
    parser.add_argument(
        '--ensembl_server',
        type=str,
        required=False,
        default='http://useast.ensembl.org/biomart',
        help='Name of the database (default: http://useast.ensembl.org/biomart)'
    )
    parser.add_argument(
        '--ref',
        type=str,
        required=False,
        default="GRCh38",
        help='The reference genome to quiery: GRCh38 (default), GRCh37, or GRCm38'
    )
    args = parser.parse_args()

    db = database.AGFusionSQlite3DB(args.name,args.ref)
    db.fetch_data(args.ensembl_server,args.p)

def main():

    parser = argparse.ArgumentParser(description='Annotate Gene Fusion (AGFusion)')
    parser.add_argument(
        '--gene5prime',
        type=str,
        required=True,
        help='5\' gene partner'
    )
    parser.add_argument(
        '--gene3prime',
        type=str,
        required=True,
        help='3\' gene partner'
    )
    parser.add_argument(
        '--junction5prime',
        type=str,
        required=True,
        help='Genomic location of predicted fuins for the 5\' gene partner'
    )
    parser.add_argument(
        '--junction3prime',
        type=str,
        required=True,
        help='Genomic location of predicted fuins for the 3\' gene partner'
    )
    parser.add_argument(
        '--db',
        type=str,
        required=True,
        help='The SQLite3 database.'
    )
    parser.add_argument(
        '--ref',
        type=str,
        required=True,
        help='Reference genome (GRCh38, GRCh37, or GRCm38)'
    )
    parser.add_argument(
        '--type',
        type=str,
        required=False,
        default='png',
        help='Image file type (png, jpeg, pdf). Default: png'
    )
    parser.add_argument(
        '--scale',
        type=int,
        required=False,
        default=-1,
        help='Length in amino acids to scale the gene fusion image (default: max length of fusion product)'
    )
    parser.add_argument(
        '--WT',
        action='store_true',
        required=False,
        help='Include this to plot wild-type architechtures of the 5\' and 3\' genes'
    )

    args = parser.parse_args()

    db = database.AGFusionSQlite3DB(args.db)

    gene5prime = model.Gene(args.gene5prime,args.junction5prime,db)

    gene3prime = model.Gene(args.gene3prime,args.junction3prime,db)

    fusion = model.Fusion(gene5prime,gene3prime)

    fusion.save_image()

    if args.WT:
        gene5prime.save_image()
        gene3prime.save_image()
