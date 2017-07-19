import sqlite3
from os.path import split, exists, join, expanduser
from os import mkdir, remove
from sys import exit
import argparse
import logging
import gzip
import shutil
from future.standard_library import install_aliases
install_aliases()
from urllib.request import urlopen
from urllib.error import HTTPError

import agfusion
from agfusion import exceptions
import pyensembl

from agfusion.utils import AGFUSION_DB_URL, AVAILABLE_ENSEMBL_SPECIES, GENOME_SHORTCUTS

def list_available_databases():
    print('\n')
    print('{:<10}\t\t{:<5}\t\t{:<20}'.format('Species','Release','Shortcut(s)'))
    for species, releases in AVAILABLE_ENSEMBL_SPECIES.items():
        for release in releases:
            shortcut = []
            for ss, data in GENOME_SHORTCUTS.items():
                if species in data and release in data:
                    shortcut.append(ss)
            print('{:<10}\t\t{:<5}\t\t{:<20}'.format(species,release,','.join(shortcut)))
    exit()

def downloaddb(args):

    if args.genome is not None:
        if args.genome not in GENOME_SHORTCUTS:
            print('Invalid genome shortcut! Use -a to see available shortcuts.')
            exit()
        else:
            species = GENOME_SHORTCUTS[args.genome][0]
            release = str(GENOME_SHORTCUTS[args.genome][1])
    else:
        if args.species is None or args.release is None:
            print("Specify --species and --release or --genome!")
            exit()
        species = args.species
        release = str(args.release)

    file_path = join(
        args.dir,
        'agfusion.' + species + '.' + release + '.db.gz')

    print("Downloading the AGFusion database to %s..." % file_path)

    db_url = AGFUSION_DB_URL + species + '.' + release + '.db.gz'

    try:
        response = urlopen(db_url)
    except HTTPError:
        print("Was unable to downloade the file %s!" % db_url)
        exit()

    fout = open(file_path,'wb')
    fout.write(response.read())
    fout.close()

    with gzip.open(file_path, 'rb') as f_in, file(file_path.replace('.gz',''), 'w') as f_out:
        shutil.copyfileobj(f_in, f_out)

    remove(file_path)

def annotate(gene5prime,junction5prime,gene3prime,junction3prime,
             outdir,colors,rename,scale,db,pyensembl_data,args):


    fusion = agfusion.Fusion(
        gene5prime=gene5prime,
        gene5primejunction=junction5prime,
        gene3prime=gene3prime,
        gene3primejunction=junction3prime,
        db=db,
        pyensembl_data=pyensembl_data,
        protein_databases=args.protein_databases,
        noncanonical=args.noncanonical
    )

    fusion.save_transcript_cdna(
        out_dir=outdir,
        middlestar=args.middlestar
    )
    fusion.save_transcript_cds(
        out_dir=outdir,
        middlestar=args.middlestar
    )
    fusion.save_proteins(
        out_dir=outdir,
        middlestar=args.middlestar
    )

    fusion.save_images(
        out_dir=outdir,
        scale=scale,
        colors=colors,
        rename=rename,
        fontsize=args.fontsize,
        height=args.height,
        width=args.width,
        dpi=args.dpi,
        no_domain_labels=args.no_domain_labels,
        plot_WT=args.WT,
        exclude=args.exclude_domain
        )
    fusion.save_tables(out_dir=outdir)

def builddb(args):

    db = agfusion.AGFusionDBBManager(
        args.dir,
        args.species,
        args.release,
        args.pfam,
        args.server
    )

    db.logger.info('Fetching alternative gene names...')

    db.fetch_gene_names()

    db.logger.info('Fetching transcript tables...')

    db.fetch_transcript_table()

    db.fetch_refseq_table()

    db.logger.info('Fetching protein annotation data...')

    db.fetch_protein_annotation()

def add_common_flags(parser):
    parser.add_argument(
        '--db',
        type=str,
        required=True,
        help='Path to the AGFusion database (e.g. --db /path/to/agfusion.homo_sapiens.87.db)'
    )
    parser.add_argument(
        '-o',
        '--out',
        type=str,
        required=True,
        help='Directory to save results'
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
        default=['pfam', 'tmhmm'],
        help='(Optional) Space-delimited list of one or more protein ' +
             'feature databases to include when visualizing proteins. ' +
             'Options are: pfam, smart, superfamily, tigrfam, pfscan (Prosite_profiles), ' +
             'tmhmm (i.e. transmembrane), seg (low_complexity regions), ncoils (coiled coil regions), prints ' +
             'pirsf, and signalp (signal peptide regions) ' +
             '(default includes pfam and tmhmm).'
    )
    parser.add_argument(
        '--recolor',
        type=str,
        required=False,
        default=None,
        action='append',
        help='(Optional) Re-color a domain. Provide the original name of the domain then your color (semi-colon delimited, all in quotes). Can specify --recolor multiples for each domain. (e.g. --color \"Pkinase_Tyr;blue\" --color \"I-set;#006600\").'
    )
    parser.add_argument(
        '--rename',
        type=str,
        required=False,
        default=None,
        action='append',
        help='(Optional) Rename a domain. Provide the original name of the domain then your new name (semi-colon delimited, all in quotes). Can specify --rename multiples for each domain. (e.g. --rename \"Pkinase_Tyr;Kinase\").'
    )
    parser.add_argument(
        '--exclude_domain',
        type=str,
        required=False,
        default=[],
        nargs='+',
        help='(Optional) Exclude a certain domain(s) from plotting by providing a space-separated list of domain names.'
    )
    parser.add_argument(
        '--type',
        type=str,
        required=False,
        default='png',
        help='(Optional) Image file type (png, jpeg, pdf). Default: png'
    )
    parser.add_argument(
        '-w',
        '--width',
        type=int,
        required=False,
        default=10,
        help='(Optional) Image width in inches (default 10).'
    )
    parser.add_argument(
        '-ht',
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
        '--debug',
        default=False,
        action='store_true',
        help='(Optional) Enable debugging logging.'
    )

def main():

    parser = argparse.ArgumentParser(
        description='Annotate Gene Fusion (AGFusion)'
    )
    subparsers = parser.add_subparsers(help='AGFusion programs.',dest="subparser_name")

    annotate_parser = subparsers.add_parser('annotate', help='Annotate and visualize a single fusion.')
    annotate_parser.add_argument(
        '-g5',
        '--gene5prime',
        type=str,
        required=True,
        help='5\' gene partner'
    )
    annotate_parser.add_argument(
        '-g3',
        '--gene3prime',
        type=str,
        required=True,
        help='3\' gene partner'
    )
    annotate_parser.add_argument(
        '-j5',
        '--junction5prime',
        type=int,
        required=True,
        help='Genomic location of predicted fuins for the 5\' gene partner. ' +
             'The 1-based position that is the last nucleotide included in ' +
             'the fusion before the junction.'
    )
    annotate_parser.add_argument(
        '-j3',
        '--junction3prime',
        type=int,
        required=True,
        help='Genomic location of predicted fuins for the 3\' gene partner. ' +
             'The 1-based position that is the first nucleotide included in ' +
             'the fusion after the junction.'
    )
    add_common_flags(annotate_parser)
    annotate_parser.add_argument(
        '--scale',
        type=int,
        required=False,
        default=-1,
        help='(Optional) Set maximum width (in amino acids) of the figure to rescale the fusion (default: max length of fusion product)'
    )

    # batch file parser

    batch_parser = subparsers.add_parser('batch', help='Annotate fusions from an output file from a fusion finding algorithm.')
    batch_parser.add_argument(
        '--file',
        type=str,
        required=True,
        help='Output file from fusion-finding algorithm.'
    )
    batch_parser.add_argument(
        '-a',
        '--algorithm',
        type=str,
        required=True,
        help='The fusion-finding algorithm. Can be one of the following: ' +
             'fusioncatcher, starfusion, or tophatfusion. (Will support more algorithms soon)'
    )
    add_common_flags(batch_parser)

    # download database

    database_parser = subparsers.add_parser('download', help='Download database for a reference genome.')
    database_parser.add_argument(
        '-d',
        '--dir',
        type=str,
        default='',
        help='(Optional) Directory to the database will be downloaded to (defaults to current working directory).'
    )
    database_parser.add_argument(
        '-g',
        '--genome',
        type=str,
        default=None,
        help='Specify the genome shortcut (e.g. hg19). To see all available ' + \
             'shortcuts run \'agfusion download -a\'. Either specify this or ' + \
             '--species and --release. '
    )
    database_parser.add_argument(
        '-s',
        '--species',
        type=str,
        default=None,
        help='The species (e.g. homo_sapiens).'
    )
    database_parser.add_argument(
        '-r',
        '--release',
        type=int,
        default=None,
        help='The ensembl release (e.g. 87).'
    )
    database_parser.add_argument(
        '-a',
        '--available',
        action='store_true',
        required=False,
        help='List available species and ensembl releases.'
    )

    # build database parser

    build_database_parser = subparsers.add_parser('build', help='Build database for a reference genome.')
    build_database_parser.add_argument(
        '--dir',
        type=str,
        required=True,
        help='Directory to write database file to.'
    )
    build_database_parser.add_argument(
        '--species',
        type=str,
        required=True,
        help='The species (e.g. homo_sapiens).'
    )
    build_database_parser.add_argument(
        '--release',
        type=int,
        required=True,
        help='The ensembl release (e.g. 87).'
    )
    build_database_parser.add_argument(
        '--pfam',
        type=str,
        required=True,
        help='File containing PFAM ID mappings.'
    )
    build_database_parser.add_argument(
        '--server',
        type=str,
        required=False,
        default='ensembldb.ensembl.org',
        help='(optional) Ensembl server (default ensembldb.ensembl.org)'
    )

    # agfusion version number

    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=agfusion.__version__
    )
    args = parser.parse_args()

    if args.subparser_name == 'build':
        builddb(args)
    elif args.subparser_name == 'download':
        if args.available:
            list_available_databases()
        else:
            downloaddb(args)
    else:
        if not exists(args.out):
            mkdir(args.out)

        # if user does not specify a sqlite database then use the one provided
        # by the package

        db_file = split(args.db)[1]
        species = db_file.split('.')[1]
        release = db_file.split('.')[2]

        assert species in AVAILABLE_ENSEMBL_SPECIES, 'unsupported species!'

        db = agfusion.AGFusionDB(args.db,debug=args.debug)
        db.build = species + '_' + str(release)

        # get the pyensembl data

        pyensembl_data = pyensembl.EnsemblRelease(release, species)

        try:
            pyensembl_data.db
        except ValueError:
            db.logger.error(
                "Missing pyensembl data. Run pyensembl install --release %s --species %s" %
                (release,species)
            )
            exit()

        colors = {}
        rename = {}

        if args.recolor is not None:
            for i in args.recolor:
                pair = i.split(';')

                assert len(pair) == 2, " did not properly specify --colors"

                if pair[0] in colors:
                    db.logger.warn("You specified colors for %s twice." % pair[0])

                colors[pair[0]] = pair[1]

        if args.rename is not None:
            for i in args.rename:
                pair = i.split(';')

                assert len(pair) == 2, " did not properly specify --rename"

                if pair[0] in rename:
                    db.logger.warn("WARNING - you rename %s twice." % pair[0])

                rename[pair[0]] = pair[1]

        if args.subparser_name == 'annotate':
            annotate(
                gene5prime=args.gene5prime,
                junction5prime=args.junction5prime,
                gene3prime=args.gene3prime,
                junction3prime=args.junction3prime,
                outdir=args.out,
                colors=colors,
                rename=rename,
                scale=args.scale,
                db=db,
                pyensembl_data=pyensembl_data,
                args=args
            )
        else:
            if args.algorithm in agfusion.parsers:
                for fusion in agfusion.parsers[args.algorithm](args.file):

                    outdir = join(
                        args.out,
                        fusion['alternative_name_5prime'] + '-' +
                        str(fusion['junction_5prime']) + '_' +
                        fusion['alternative_name_3prime'] + '-' +
                        str(fusion['junction_3prime'])
                    )

                    if not exists(outdir):
                        mkdir(outdir)
                    else:
                        db.logger.warn('The following output directory already exists! %s' % outdir)

                    try:
                        annotate(
                            gene5prime=fusion['ensembl_5prime'],
                            junction5prime=fusion['junction_5prime'],
                            gene3prime=fusion['ensembl_3prime'],
                            junction3prime=fusion['junction_3prime'],
                            outdir=outdir,
                            colors=colors,
                            rename=rename,
                            scale=None,
                            db=db,
                            pyensembl_data=pyensembl_data,
                            args=args
                        )
                    except exceptions.GeneIDException3prime:
                        db.logger.warn("No Ensembl ID found for {0}! Check its spelling and if you are using the right genome build.".format(fusion['ensembl_5prime']))
                    except exceptions.GeneIDException3prime:
                        db.logger.warn("No Ensembl ID found for {0}! Check its spelling and if you are using the right genome build.".format(fusion['ensembl_3prime']))
                    except exceptions.JunctionException5prime:
                        db.logger.warn("No Ensembl ID found for {0}! Check its spelling and if you are using the right genome build.".format(fusion['ensembl_5prime']))
                    except exceptions.JunctionException3prime:
                        db.logger.warn("No Ensembl ID found for {0}! Check its spelling and if you are using the right genome build.".format(fusion['ensembl_3prime']))
                    except exceptions.TooManyGenesException as e:
                        db.logger.warn("Multiple Ensembl IDs found matching one of your genes.")
