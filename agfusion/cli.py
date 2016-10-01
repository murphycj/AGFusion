import sqlite3
import os
import sys
import argparse
import logging
import urllib

import agfusion

def build_db():

    parser = argparse.ArgumentParser(
        description='Build the SQLite3 database for a reference ' + \
        'genomes by querying Biomart. The the database given by --database ' + \
        'already exists then that portion will be overwritten.'
    )
    parser.add_argument(
        '--database',
        type=str,
        required=True,
        help='Path to the database file (e.g. agfusion.db)'
    )
    parser.add_argument(
        '--genome',
        type=str,
        required=True,
        help='GRCh38, GRCh37, or GRCm38'
    )
    parser.add_argument(
        '--p',
        type=int,
        default=1,
        help='(Optional) Number of processers to use to fetch data from Biomart (default 1).'
    )

    args = parser.parse_args()

    db = agfusion.AGFusionDBBManager(args.database,args.genome)

    db.logger.info('Fetching data from Biomart...')

    db.fetch_data(args.p)

    db.logger.info('Retrieving PFAM domain name mapping file...')

    pfam_file = os.path.join(os.path.dirname(args.database),'pdb_pfam_mapping.txt')
    r = urllib.urlretrieve('ftp://ftp.ebi.ac.uk/pub/databases/Pfam//mappings/pdb_pfam_mapping.txt',pfam_file)

    db.logger.info('Adding PFAM data to the database...')

    db.add_pfam(pfam_file=pfam_file)

    db.logger.info('Done!')

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
        '--genome',
        type=str,
        required=True,
        help='GRCh38, GRCh37, or GRCm38'
    )
    parser.add_argument(
        '--colors',
        type=str,
        required=False,
        nargs='+',
        default=None,
        help='(Optional) Space-delimited list of domain name and color to ' + \
            'specify certain colors for domains. Format --color domain_name:color' + \
            ' (e.g. --color Pkinase_Tyr:blue I-set:#006600). ' + \
            'Can use specific color names for hex representation. Default ' + \
            'blue for everything.'
    )
    parser.add_argument(
        '--rename',
        type=str,
        required=False,
        nargs='+',
        default=None,
        help='(Optional) Space-delimited list of domain name and new name to ' + \
            'rename particular domains. Format --rename domain_name:new_domain_name' + \
            ' (e.g. --rename Pkinase_Tyr:Kinase).'
    )
    parser.add_argument(
        '--type',
        type=str,
        required=False,
        default='png',
        help='(Optional) Image file type (png, jpeg, pdf). Default: png'
    )
    parser.add_argument(
        '--scale',
        type=int,
        required=False,
        default=-1,
        help='(Optional) Length in amino acids to scale the gene fusion image (default: max length of fusion product)'
    )
    parser.add_argument(
        '--WT',
        action='store_true',
        required=False,
        help='(Optional) Include this to plot wild-type architechtures of the 5\' and 3\' genes'
    )
    parser.add_argument(
        '--middlestar',
        action='store_true',
        required=False,
        help='(Optional) Insert a * at the junction position for the cdna, cds, and protein sequences (default False).'
    )

    args = parser.parse_args()

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    db = agfusion.AGFusionDB(args.db,args.genome)

    gene5prime = agfusion.Gene(
        gene=args.gene5prime,
        junction=args.junction5prime,
        db=db
    )

    gene3prime = agfusion.Gene(
        gene=args.gene3prime,
        junction=args.junction3prime,
        db=db
    )

    fusion = agfusion.Fusion(
        gene5prime=gene5prime,
        gene3prime=gene3prime,
        db=db,
        middlestar=args.middlestar
    )

    fusion.save_transcript_cdna(args.out)
    fusion.save_transcript_cds(args.out)
    fusion.save_proteins(args.out)

    colors={}
    rename = {}

    if args.colors is not None:
        for i in args.colors:
            pair = i.split(':')

            assert len(pair)==2," did not properly specify --colors"

            if pair[0] in colors:
                print "!!! WARNING - you specified colors for %s twice." % pair[0]

            colors[pair[0]] = pair[1]

    if args.rename is not None:
        for i in args.rename:
            pair = i.split(':')

            assert len(pair)==2," did not properly specify --rename"

            if pair[0] in rename:
                print "!!! WARNING - you rename %s twice." % pair[0]

            rename[pair[0]] = pair[1]

    fusion.save_images(out_dir=args.out,colors=colors,rename=rename)

    if args.WT:
        gene5prime.save_image(args.out)
        gene3prime.save_image(args.out)
