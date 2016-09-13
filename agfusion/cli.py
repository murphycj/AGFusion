import sqlite3
import os
import sys
import argparse

import agfusion
import pyensembl

def build_db():

    parser = argparse.ArgumentParser(
        description='Build or update the SQLite3 database for a reference ' + \
        'genomes by querying Biomart. If the database given by --name already ' + \
        'exists, then species-specific portion of the database will be updated.'
    )
    parser.add_argument(
        '--database',
        type=str,
        required=True,
        help='Path to the database file (e.g. agfusion.db)'
    )
    parser.add_argument(
        '--ensembl_dataset',
        type=str,
        required=True,
        help='The reference ensembl dataset to query biomart (hsapiens_gene_ensembl or mmusculus_gene_ensembl)'
    )
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        help='Human or mouse'
    )
    parser.add_argument(
        '--release',
        type=str,
        required=True,
        help='Ensembl release'
    )
    parser.add_argument(
        '--ensembl_server',
        type=str,
        required=False,
        default='http://useast.ensembl.org/biomart',
        help='(Optional) Name of the database (default: http://useast.ensembl.org/biomart)'
    )
    parser.add_argument(
        '--p',
        type=int,
        default=1,
        help='(Optional) Number of processers to use to fetch data from Biomart (default 1).'
    )
    args = parser.parse_args()

    data = pyensembl.EnsemblRelease(args.release,args.species)

    if args.species.lower()=='human':
        dataset='hsapiens_gene_ensembl'
    elif args.species.lower()=='mouse':
        dataset='mmusculus_gene_ensembl'

    db = agfusion.AGFusionDBBManager(args.database)
    db.fetch_data(args.ensembl_server,dataset,args.p,data.transcript_ids())
    db.add_pfam()

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
        type=int,
        required=True,
        help='Genomic location of predicted fuins for the 5\' gene partner'
    )
    parser.add_argument(
        '--junction3prime',
        type=int,
        required=True,
        help='Genomic location of predicted fuins for the 3\' gene partner'
    )
    parser.add_argument(
        '--out',
        type=str,
        required=True,
        help='Directory to save results'
    )
    parser.add_argument(
        '--db',
        type=str,
        required=True,
        help='The SQLite3 database.'
    )
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        help='Reference genome (GRCh38, GRCh37, or GRCm38)'
    )
    parser.add_argument(
        '--release',
        type=str,
        required=True,
        help='Ensembl release'
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

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    db = agfusion.AGFusionDB(args.db)
    data = pyensembl.EnsemblRelease(args.release,args.species)

    gene5prime = agfusion.Gene(
        gene=data.gene_by_id(args.gene5prime),
        junction=args.junction5prime,
        db=db
    )

    gene3prime = agfusion.Gene(
        gene=data.gene_by_id(args.gene3prime),
        junction=args.junction3prime,
        db=db
    )

    fusion = agfusion.Fusion(gene5prime,gene3prime,db=db)

    fusion.save_transcript_cdna(args.out)
    fusion.save_transcript_cds(args.out)
    fusion.save_proteins(args.out)
    fusion.save_image(args.out)

    if args.WT:
        gene5prime.save_image(args.out)
        gene3prime.save_image(args.out)
