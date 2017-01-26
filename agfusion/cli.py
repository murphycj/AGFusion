import sqlite3
import os
import sys
import argparse
import logging
import urllib

import agfusion
import pyensembl

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
        help='Genomic location of predicted fuins for the 5\' gene partner. ' + \
             'The 1-based position that is the last nucleotide included in ' + \
             'the fusion before the junction.'
    )
    parser.add_argument(
        '--junction3prime',
        type=int,
        required=True,
        help='Genomic location of predicted fuins for the 3\' gene partner. ' + \
             'The 1-based position that is the first nucleotide included in ' + \
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
        help='GRCh38, GRCh37, or GRCm38'
    )
    parser.add_argument(
        '--db',
        type=str,
        default=None,
        required=False,
        help='(Optional) The SQLite3 database. Defaults to using the database provided by the package.'
    )
    parser.add_argument(
        '--protein_databases',
        type=str,
        required=False,
        nargs='+',
        default=['pfam','tmhmm'],
        help='(Optional) Space-delimited list of one or more protein ' + \
             'feature databases to include when visualizing proteins. ' + \
             'Options are: pfam, tigrfam, prints, hmmpanther, ' + \
             'blastprodom, gene3d, hamap, pirsf, ncoils, superfamily, seg, ' + \
             'signalp, scanprosite, pfscan, tmhmm, and smart. ' + \
             '(default includes pfam and tmhmm).'
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
    #parser.add_argument(
    #    '--WT',
    #    action='store_true',
    #    required=False,
    #    help='(Optional) Include this to plot wild-type architechtures of the 5\' and 3\' genes'
    #)
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

    #if user does not specify a sqlite database then use the one provided
    #by the package

    if args.db is None:
        file_path = os.path.split(__file__)[0]
        db = agfusion.AGFusionDB(os.path.join(file_path,'data','agfusion.db'))
    else:
        db = agfusion.AGFusionDB(args.db)

    #get the pyensembl data

    if args.genome=='GRCm38':
        pyensembl_data = pyensembl.EnsemblRelease(84,'mouse')
    elif args.genome=='GRCh38':
        pyensembl_data = pyensembl.EnsemblRelease(84,'human')
    elif args.genome=='GRCh37':
        pyensembl_data = pyensembl.EnsemblRelease(75,'human')
    else:
        db.logger.error(' You provided an incorrect reference genome. Use one of the following: GRCh38, GRCh37, or GRCm38')
        sys.exit()

    fusion = agfusion.Fusion(
        gene5prime=args.gene5prime,
        gene5primejunction=args.junction5prime,
        gene3prime=args.gene3prime,
        gene3primejunction=args.junction3prime,
        db=db,
        pyensembl_data=pyensembl_data,
        protein_databases=args.protein_databases
    )

    fusion.save_transcript_cdna(out_dir=args.out,middlestar=args.middlestar)
    fusion.save_transcript_cds(out_dir=args.out,middlestar=args.middlestar)
    fusion.save_proteins(out_dir=args.out,middlestar=args.middlestar)

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

    fusion.save_images(
        out_dir=args.out,
        scale=args.scale,
        colors=colors,
        rename=rename,
        fontsize=args.fontsize,
        height=args.height,
        width=args.width,
        dpi=args.dpi,
        no_domain_labels=args.no_domain_labels
        )

#    if args.WT:
#        gene5prime.save_image(args.out)
#        gene3prime.save_image(args.out)
