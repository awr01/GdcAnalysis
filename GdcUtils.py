import argparse
from GdcUtils import *

parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')

parser.add_argument( '--genes' , nargs='+' , default=[] , help='Genes of interest')
parser.add_argument( '--mutations' , nargs='+' , default=[] , help='Mutations of interest')

parser.add_argument( '--gene-file' , help='Genes of interest')
parser.add_argument( '--mutation-file' , help='Mutations of interest')

args = parser.parse_args()

if not args.src.endswith( ".tgz" ): raise Exception( "Source file must have '.tgz' file-extension" )
if not ( args.dest.endswith( ".tsv" ) or args.dest.endswith( ".xlsx" ) ): raise Exception( "Destination file must have '.tsv' or '.xlsx' file-extension" )

if not args.gene_file is None:
  with open( args.gene_file , "r" ) as src:
    for line in src: args.genes.append( line.strip() )

if not args.mutation_file is None:
  with open( args.mutation_file , "r" ) as src:
    for line in src: args.mutations.append( line.strip() )
