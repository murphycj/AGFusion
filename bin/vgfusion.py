import argparse

from vgfusion import *

parser = argparse.ArgumentParser(description='Visualize Gene Fusion (VGFusion)')
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
    '--5primeLocation',
    type=str,
    required=True,
    help='Genomic location of predicted fuins for the 5\' gene partner'
)
parser.add_argument(
    '--3primeLocation',
    type=str,
    required=True,
    help='Genomic location of predicted fuins for the 3\' gene partner'
)
parser.add_argument(
    '--ref',
    type=str,
    required=True,
    help='Reference genome (GRCh38, GRCh37, or GRCm38)'
)
parser.add_argument(
    '--database',
    type=str,
    required=False,
    default=None,
    help='Path to database. If not specified VGFusion will make RESTful API calls to Ensembl instead'
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

main(args=args)
