""" cli.py
Command line interface
"""

import argparse
import gzip
import shutil
import sys
from os import mkdir, remove
from os.path import exists, join, split
from pathlib import Path

from future.standard_library import install_aliases

install_aliases()
from urllib.error import HTTPError
from urllib.request import urlopen

import pyensembl

from agfusion import __version__, database, exceptions, model, parsers, utils


def list_available_databases():
    """
    List the available databases that can be downloaded.
    #"""

    print("\n")
    print("Species\t\tRelease\t\tShortbut(s)")
    for species, releases in utils.AVAILABLE_ENSEMBL_SPECIES.items():
        for release in releases:
            shortcut = []
            for genome, data in utils.GENOME_SHORTCUTS.items():
                if species in data and release in data:
                    shortcut.append(genome)
            shortcut = ",".join(shortcut)
            print("{species:<10}\t\t{release:<5}\t\t{shortcut:<20}")
    sys.exit(0)


def downloaddb(args):
    """
    Download the AGFusion database from github
    """

    if args.genome is not None:
        if args.genome not in utils.GENOME_SHORTCUTS:
            print("Invalid genome shortcut! Use -a to see available shortcuts.")
            sys.exit(1)
        else:
            species = utils.GENOME_SHORTCUTS[args.genome][0]
            release = str(utils.GENOME_SHORTCUTS[args.genome][1])
    else:
        if args.species is None or args.release is None:
            print("Specify --species and --release or --genome!")
            sys.exit(1)
        species = args.species
        release = str(args.release)

    file_path = join(args.dir, "agfusion." + species + "." + release + ".db.gz")

    print(f"Downloading the AGFusion database to {file_path}...")

    db_url = utils.AGFUSION_DB_URL + species + "." + release + ".db.gz"

    try:
        response = urlopen(db_url)
    except HTTPError:
        print(f"Was unable to download the file {db_url}!")
        sys.exit(1)

    fout = open(file_path, "wb")
    fout.write(response.read())
    fout.close()

    with gzip.open(file_path, "rb") as f_in, open(file_path.replace(".gz", ""), "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    remove(file_path)


def annotate(
    gene5prime,
    junction5prime,
    gene3prime,
    junction3prime,
    agfusion_db,
    pyensembl_data,
    args,
    outdir=None,
    colors=None,
    rename=None,
    scale=None,
    batch_out_dir=None,
):
    """
    Annotate the gene fusion
    """

    fusion = model.Fusion(
        gene5prime=gene5prime,
        gene5primejunction=junction5prime,
        gene3prime=gene3prime,
        gene3primejunction=junction3prime,
        db=agfusion_db,
        pyensembl_data=pyensembl_data,
        protein_databases=args.protein_databases,
        noncanonical=args.noncanonical,
    )

    if batch_out_dir is not None:

        outdir = join(
            batch_out_dir,
            fusion.gene5prime.gene.name
            + "-"
            + str(junction5prime)
            + "_"
            + fusion.gene3prime.gene.name
            + "-"
            + str(junction3prime),
        )

    fusion.save_transcript_cdna(out_dir=outdir, middlestar=args.middlestar)
    fusion.save_transcript_cds(out_dir=outdir, middlestar=args.middlestar)
    fusion.save_proteins(out_dir=outdir, middlestar=args.middlestar)

    fusion.save_images(
        out_dir=outdir,
        file_type=args.type,
        scale=scale,
        colors=colors,
        rename=rename,
        fontsize=args.fontsize,
        height=args.height,
        width=args.width,
        dpi=args.dpi,
        no_domain_labels=args.no_domain_labels,
        plot_wt=args.WT,
        exclude=args.exclude_domain,
    )
    fusion.save_tables(out_dir=outdir)


def batch_mode(args, agfusion_db, pyensembl_data, rename, colors):
    """
    Batch mode for annotation fusions from output from a fusion-finding
    algorithm
    """

    if not exists(args.out):
        mkdir(args.out)
    else:
        agfusion_db.logger.warn(f"Output directory {args.out} already exists! Overwriting...")

    if not Path(args.file).exists():
        FileNotFoundError(f"File not found {args.file}")
        sys.exit(1)

    if args.algorithm in parsers.parsers:
        for fusion in parsers.parsers[args.algorithm](args.file, agfusion_db.logger):

            try:
                annotate(
                    gene5prime=fusion["gene5prime"],
                    junction5prime=fusion["gene5prime_junction"],
                    gene3prime=fusion["gene3prime"],
                    junction3prime=fusion["gene3prime_junction"],
                    agfusion_db=agfusion_db,
                    pyensembl_data=pyensembl_data,
                    args=args,
                    colors=colors,
                    rename=rename,
                    scale=None,
                    batch_out_dir=args.out,
                )
            except exceptions.GeneIDException as error:
                agfusion_db.logger.error(error)
            except exceptions.JunctionException as error:
                agfusion_db.logger.error(error)
            except exceptions.TooManyGenesException as error:
                agfusion_db.logger.error(error)
    else:
        available = ",".join(parsers.parsers.keys())
        agfusion_db.logger.error(
            (
                f"'{args.algorithm}' is not an available option for -a! Choose one of the "
                + f"following: {available}."
            )
        )
        sys.exit(1)


def builddb(args):
    """
    Build a AGFusion database
    """

    agfusion_db = database.AGFusionDBBManager(
        args.dir, args.species, args.release, args.pfam, args.server
    )

    agfusion_db.logger.info("Fetching alternative gene names...")

    agfusion_db.fetch_gene_names()

    agfusion_db.logger.info("Fetching transcript tables...")

    agfusion_db.fetch_transcript_table()

    agfusion_db.fetch_refseq_table()

    agfusion_db.logger.info("Fetching protein annotation data...")

    agfusion_db.fetch_protein_annotation()


def add_common_flags(parser):
    """
    Add commaond line flags that are common to multiple sub parsers
    """

    parser.add_argument(
        "-db",
        "--database",
        type=str,
        required=True,
        help="Path to the AGFusion database (e.g. --db /path/to/agfusion.homo_sapiens.87.db)",
    )
    parser.add_argument("-o", "--out", type=str, required=True, help="Directory to save results")
    parser.add_argument(
        "-nc",
        "--noncanonical",
        action="store_true",
        required=False,
        default=False,
        help="(Optional) Include non-canonical gene transcripts "
        + "in the analysis (default False).",
    )
    parser.add_argument(
        "--protein_databases",
        type=str,
        required=False,
        nargs="+",
        default=["pfam", "tmhmm"],
        help="(Optional) Space-delimited list of one or more protein "
        + "feature databases to include when visualizing proteins. "
        + "Options are: pfam, smart, superfamily, tigrfam, pfscan (Prosite_profiles), "
        + "tmhmm (i.e. transmembrane), seg (low_complexity regions), ncoils "
        + "(coiled coil regions), prints, "
        + "pirsf, and signalp (signal peptide regions) "
        + "(default: --protein_databases pfam and tmhmm).",
    )
    parser.add_argument(
        "--recolor",
        type=str,
        required=False,
        default=None,
        action="append",
        help="(Optional) Re-color a domain. Provide the original name of "
        + "the domain then your color (semi-colon delimited, all in "
        + "quotes). Can specify --recolor multiples for each domain. "
        + '(e.g. --color "Pkinase_Tyr;blue" --color "I-set;#006600").',
    )
    parser.add_argument(
        "--rename",
        type=str,
        required=False,
        default=None,
        action="append",
        help="(Optional) Rename a domain. Provide the original name of "
        + "the domain then your new name (semi-colon delimited, "
        + "all in quotes). Can specify --rename multiples for each "
        + 'domain. (e.g. --rename "Pkinase_Tyr;Kinase").',
    )
    parser.add_argument(
        "--exclude_domain",
        type=str,
        required=False,
        default=[],
        nargs="+",
        help="(Optional) Exclude a certain domain(s) from plotting "
        + "by providing a space-separated list of domain names.",
    )
    parser.add_argument(
        "--type",
        type=str,
        required=False,
        default="png",
        help="(Optional) Image file type (png, jpeg, pdf). Default: png",
    )
    parser.add_argument(
        "-w",
        "--width",
        type=int,
        required=False,
        default=10,
        help="(Optional) Image width in inches (default 10).",
    )
    parser.add_argument(
        "-ht",
        "--height",
        type=int,
        required=False,
        default=3,
        help="(Optional) Image file height in inches (default 3).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        required=False,
        default=None,
        help="(Optional) Dots per inch.",
    )
    parser.add_argument(
        "--fontsize",
        type=int,
        required=False,
        default=12,
        help="(Optional) Fontsize (default 12).",
    )
    parser.add_argument(
        "--WT",
        action="store_true",
        required=False,
        help="(Optional) Include this to plot wild-type architechtures " + "of the 5' and 3' genes",
    )
    parser.add_argument(
        "-ms",
        "--middlestar",
        action="store_true",
        required=False,
        help="(Optional) Insert a * at the junction position for the "
        + "cdna, cds, and protein sequences (default False).",
    )
    parser.add_argument(
        "-ndl",
        "--no_domain_labels",
        action="store_true",
        required=False,
        help="(Optional) Do not label domains.",
    )
    parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="(Optional) Enable debugging logging.",
    )


def build_cli():
    """Build the CLI.

    Returns:
        argparse.Arguments: Parsed arguments.
    """

    parser = argparse.ArgumentParser(description="Annotate Gene Fusion (AGFusion)")
    subparsers = parser.add_subparsers(help="AGFusion programs.", dest="subparser_name")

    annotate_parser = subparsers.add_parser(
        "annotate", help="Annotate and visualize a single fusion."
    )
    annotate_parser.add_argument(
        "-g5", "--gene5prime", type=str, required=True, help="5' gene partner"
    )
    annotate_parser.add_argument(
        "-g3", "--gene3prime", type=str, required=True, help="3' gene partner"
    )
    annotate_parser.add_argument(
        "-j5",
        "--junction5prime",
        type=int,
        required=True,
        help="Genomic location of predicted fuins for the 5' gene partner. "
        + "The 1-based position that is the last nucleotide included in "
        + "the fusion before the junction.",
    )
    annotate_parser.add_argument(
        "-j3",
        "--junction3prime",
        type=int,
        required=True,
        help="Genomic location of predicted fuins for the 3' gene partner. "
        + "The 1-based position that is the first nucleotide included in "
        + "the fusion after the junction.",
    )
    add_common_flags(annotate_parser)
    annotate_parser.add_argument(
        "--scale",
        type=int,
        required=False,
        default=-1,
        help="(Optional) Set maximum width (in amino acids) of the "
        + "figure to rescale the fusion (default: max length of "
        + "fusion product)",
    )

    # batch file parser

    batch_parser = subparsers.add_parser(
        "batch",
        help="Annotate fusions from an output file from a fusion " + "finding algorithm.",
    )
    batch_parser.add_argument(
        "-f",
        "--file",
        type=str,
        required=True,
        help="Output file from fusion-finding algorithm.",
    )
    batch_parser.add_argument(
        "-a",
        "--algorithm",
        type=str,
        required=True,
        help="The fusion-finding algorithm. Can be one of the following: "
        + ", ".join(parsers.parsers.keys())
        + ".",
    )
    add_common_flags(batch_parser)

    # download database

    database_parser = subparsers.add_parser(
        "download", help="Download database for a reference genome."
    )
    database_parser.add_argument(
        "-d",
        "--dir",
        type=str,
        default="",
        help="(Optional) Directory to the database will be downloaded "
        + "to (defaults to current working directory).",
    )
    database_parser.add_argument(
        "-g",
        "--genome",
        type=str,
        default=None,
        help="Specify the genome shortcut (e.g. hg19). To see all"
        + "available shortcuts run 'agfusion download -a'. Either "
        + "specify this or --species and --release.",
    )
    database_parser.add_argument(
        "-s",
        "--species",
        type=str,
        default=None,
        help="The species (e.g. homo_sapiens).",
    )
    database_parser.add_argument(
        "-r", "--release", type=int, default=None, help="The ensembl release (e.g. 87)."
    )
    database_parser.add_argument(
        "-a",
        "--available",
        action="store_true",
        required=False,
        help="List available species and ensembl releases.",
    )

    # build database parser

    build_database_parser = subparsers.add_parser(
        "build", help="Build database for a reference genome."
    )
    build_database_parser.add_argument(
        "-d",
        "--dir",
        type=str,
        required=True,
        help="Directory to write database file to.",
    )
    build_database_parser.add_argument(
        "-s",
        "--species",
        type=str,
        required=True,
        help="The species (e.g. homo_sapiens).",
    )
    build_database_parser.add_argument(
        "-r",
        "--release",
        type=int,
        required=True,
        help="The ensembl release (e.g. 87).",
    )
    build_database_parser.add_argument(
        "--pfam", type=str, required=True, help="File containing PFAM ID mappings."
    )
    build_database_parser.add_argument(
        "--server",
        type=str,
        required=False,
        default="ensembldb.ensembl.org",
        help="(optional) Ensembl server (default ensembldb.ensembl.org)",
    )

    # agfusion version number

    parser.add_argument("-v", "--version", action="version", version=__version__)
    args = parser.parse_args()

    if args.subparser_name is None:
        parser.print_help()
        sys.exit(1)

    return args


def parse_names_and_colors(args, agfusion_db):
    """parse the re-coloring and re-naming

    Args:
        args (argparse.Namespace): arguments
        agfusion_db (object): AGFusion datapase pointer

    Returns:
        list[dict]: re-mapping of colors and name
    """

    colors = {}
    rename = {}

    if args.rename is not None:
        for i in args.rename:
            pair = i.split(";")

            assert len(pair) == 2, " did not properly specify --rename"

            if pair[0] in rename:
                agfusion_db.logger.warn(f"WARNING - you rename {pair[0]} twice.")

            rename[pair[0]] = pair[1]

    if args.recolor is not None:
        for i in args.recolor:
            pair = i.split(";")

            assert len(pair) == 2, " did not properly specify --colors"

            if pair[0] in colors:
                agfusion_db.logger.warn(f"You specified colors for {pair[0]} twice.")

            if pair[0] in rename:
                colors[rename[pair[0]]] = pair[1]
            else:
                colors[pair[0]] = pair[1]

    return colors, rename


def main():
    """
    Main function for processing command line options
    """

    args = build_cli()

    if args.subparser_name == "build":
        builddb(args)
        sys.exit(0)
    elif args.subparser_name == "download":
        if args.available:
            list_available_databases()
        else:
            downloaddb(args)
        sys.exit(0)

    # single or batch mode

    if not exists(args.out):
        mkdir(args.out)

    # if user does not specify a sqlite database then use the one provided
    # by the package

    db_file = split(args.database)[1]
    species = db_file.split(".")[1]
    release = db_file.split(".")[2]

    assert species in utils.AVAILABLE_ENSEMBL_SPECIES, "unsupported species!"

    agfusion_db = database.AGFusionDB(args.database, debug=args.debug)
    agfusion_db.build = species + "_" + str(release)

    # get the pyensembl data

    pyensembl_data = pyensembl.EnsemblRelease(release, species)

    try:
        pyensembl_data.db
    except ValueError:
        agfusion_db.logger.error(
            "Missing pyensembl data. Run pyensembl install --release "
            + f"{release} --species {species}"
        )
        sys.exit(1)

    colors, rename = parse_names_and_colors(args, agfusion_db)

    # check image file type is valid

    if args.type not in ["png", "pdf", "jpeg"]:
        agfusion_db.logger.error(f"ERROR - provided an incorrect image file type: {args.type}.")
        sys.exit(1)

    if args.subparser_name == "annotate":
        annotate(
            gene5prime=args.gene5prime,
            junction5prime=args.junction5prime,
            gene3prime=args.gene3prime,
            junction3prime=args.junction3prime,
            agfusion_db=agfusion_db,
            pyensembl_data=pyensembl_data,
            args=args,
            outdir=args.out,
            colors=colors,
            rename=rename,
            scale=args.scale,
        )
    elif args.subparser_name == "batch":
        batch_mode(args, agfusion_db, pyensembl_data, rename, colors)
