import sqlite3
import os
import sys
import argparse
import logging
import urllib
import gzip
import shutil

import agfusion
import pyensembl


def builddb():

    parser = argparse.ArgumentParser(
        description='Build the SQLite3 database for a reference ' +
        'genomes by querying Biomart. The the database given by --database ' +
        'already exists then that portion will be overwritten.'
    )
    parser.add_argument(
        '--database',
        type=str,
        required=True,
        help='Path to the database file (e.g. agfusion.db)'
    )
    parser.add_argument(
        '--build',
        type=str,
        required=True,
        help='homo_sapiens_core_84_38 (for GRCh38), ' +
        'homo_sapiens_core_75_37 (for GRCh37), or ' +
        'mus_musculus_core_84_38 (for GRCm38)'
    )
    parser.add_argument(
        '--server',
        type=str,
        required=False,
        default='ensembldb.ensembl.org',
        help='(optional) Ensembl server (default ensembldb.ensembl.org)'
    )
    args = parser.parse_args()

    db = agfusion.AGFusionDBBManager(args.database, args.build, args.server)

    db.logger.info('Fetching alternative gene names...')

    db.fetch_gene_names()

    db.logger.info('Fetching transcript tables...')

    db.fetch_transcript_table()

    db.logger.info('Fetching protein annotation data...')

    db.fetch_protein_annotation()


def main():

    parser = argparse.ArgumentParser(
        description='Annotate Gene Fusion (AGFusion)'
    )
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
        help='Genomic location of predicted fuins for the 5\' gene partner. ' +
             'The 1-based position that is the last nucleotide included in ' +
             'the fusion before the junction.'
    )
    parser.add_argument(
        '--junction3prime',
        type=int,
        required=True,
        help='Genomic location of predicted fuins for the 3\' gene partner. ' +
             'The 1-based position that is the first nucleotide included in ' +
             'the fusion after the junction.'
    )
    parser.add_argument(
        '--out',
        type=str,
        required=True,
        help='Directory to save results'
    )
    parser.add_argument(
        '--genome',
        type=str,
        required=True,
        help='GRCh38 (or homo_sapiens_core_84_38), GRCh37 (or homo_sapiens_core_75_37), or GRCm38 (or mus_musculus_core_84_38)'
    )
    parser.add_argument(
        '--db',
        type=str,
        default=None,
        required=False,
        help='(Optional) The SQLite3 database. Defaults to using the ' +
             'database provided by the package.'
    )
    parser.add_argument(
        '--noncanonical',
        action='store_true',
        required=False,
        default=False,
        help='(Optional) Include non-canonical gene transcripts ' +
             'in the analysis (default False).'
    )
    parser.add_argument(
        '--protein_databases',
        type=str,
        required=False,
        nargs='+',
        default=['Pfam', 'transmembrane'],
        help='(Optional) Space-delimited list of one or more protein ' +
             'feature databases to include when visualizing proteins. ' +
             'Options are: Pfam, Smart, Superfamily, TIGRfam, Prosite_profiles, ' +
             'transmembrane, low_complexity, coiled_coil, Prints ' +
             'PIRSF, and signal_peptide ' +
             '(default includes Pfam and transmembrane).'
    )
    parser.add_argument(
        '--colors',
        type=str,
        required=False,
        default=None,
        action='append',
        help='(Optional) Re-color a domain. Provide the original name of the domain then your color (semi-colon delimited, all in quotes). Can specify --recolor multiples for each domain. (e.g. --color \"Pkinase_Tyr:blue\" --color \"I-set:#006600\").'
    )
    parser.add_argument(
        '--rename',
        type=str,
        required=False,
        default=None,
        action='append',
        help='(Optional) Rename a domain. Provide the original name of the domain then your new name (semi-colon delimited, all in quotes). Can specify --rename multiples for each domain. (e.g. --rename \"Pkinase_Tyr:Kinase\").'
    )
    parser.add_argument(
        '--type',
        type=str,
        required=False,
        default='png',
        help='(Optional) Image file type (png, jpeg, pdf). Default: png'
    )
    parser.add_argument(
        '--width',
        type=int,
        required=False,
        default=10,
        help='(Optional) Image width in inches (default 10).'
    )
    parser.add_argument(
        '--height',
        type=int,
        required=False,
        default=3,
        help='(Optional) Image file height in inches (default 3).'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        required=False,
        default=None,
        help='(Optional) Dots per inch.'
    )
    parser.add_argument(
        '--fontsize',
        type=int,
        required=False,
        default=12,
        help='(Optional) Fontsize (default 12).'
    )
    parser.add_argument(
        '--scale',
        type=int,
        required=False,
        default=None,
        help='(Optional) Set maximum width (in amino acids) of the figure to rescale the fusion (default: max length of fusion product)'
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
    parser.add_argument(
        '--no_domain_labels',
        action='store_true',
        required=False,
        help='(Optional) Do not label domains.'
    )
    parser.add_argument(
        '--version',
        action='version',
        version=agfusion.__version__
    )

    args = parser.parse_args()

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    # if user does not specify a sqlite database then use the one provided
    # by the package

    if args.db is None:
        file_path = os.path.join(
            os.path.split(__file__)[0],
            'data',
            'agfusion.db'
        )

        db = agfusion.AGFusionDB(
            file_path
        )
    else:
        db = agfusion.AGFusionDB(args.db)

    # get the pyensembl data

    if args.genome == 'GRCm38' or args.genome == 'mus_musculus_core_84_38':
        pyensembl_data = pyensembl.EnsemblRelease(84, 'mouse')
        db.build = 'mus_musculus_core_84_38'
    elif args.genome == 'GRCh38' or args.genome == 'homo_sapiens_core_84_38':
        pyensembl_data = pyensembl.EnsemblRelease(84, 'human')
        db.build = 'homo_sapiens_core_84_38'
    elif args.genome == 'GRCh37' or args.genome == 'homo_sapiens_core_75_37':
        pyensembl_data = pyensembl.EnsemblRelease(75, 'human')
        db.build = 'homo_sapiens_core_75_37'
    else:
        db.logger.error(
            'You provided an incorrect reference genome. '
            'Use one of the following: GRCh38, GRCh37, or GRCm38'
        )
        sys.exit()

    fusion = agfusion.Fusion(
        gene5prime=args.gene5prime,
        gene5primejunction=args.junction5prime,
        gene3prime=args.gene3prime,
        gene3primejunction=args.junction3prime,
        db=db,
        pyensembl_data=pyensembl_data,
        protein_databases=args.protein_databases,
        noncanonical=args.noncanonical
    )

    fusion.save_transcript_cdna(
        out_dir=args.out,
        middlestar=args.middlestar
    )
    fusion.save_transcript_cds(
        out_dir=args.out,
        middlestar=args.middlestar
    )
    fusion.save_proteins(
        out_dir=args.out,
        middlestar=args.middlestar
    )

    colors = {}
    rename = {}

    if args.colors is not None:
        for i in args.colors:
            pair = i.split(';')

            assert len(pair) == 2, " did not properly specify --colors"

            if pair[0] in colors:
                print "!!! WARNING - you specified colors for %s twice." % pair[0]

            colors[pair[0]] = pair[1]

    if args.rename is not None:
        for i in args.rename:
            pair = i.split(';')

            assert len(pair) == 2, " did not properly specify --rename"

            if pair[0] in rename:
                print "!!! WARNING - you rename %s twice." % pair[0]

            rename[pair[0]] = pair[1]

    fusion.save_images(
        out_dir=args.out,
        scale=args.scale,
        colors=colors,
        rename=rename,
        fontsize=args.fontsize,
        height=args.height,
        width=args.width,
        dpi=args.dpi,
        no_domain_labels=args.no_domain_labels,
        plot_WT=args.WT
        )
