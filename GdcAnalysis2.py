import argparse , requests, json, tarfile, io, codecs, gzip , tqdm , os , hashlib , bz2 , _pickle

from GdcLib2 import *

# Read the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--prefix' , required=True , help='The prefix that will be used for all files')
args = parser.parse_args()


lCases = LoadCases( f"{args.prefix}.tgz" )
  
  
for i,j in lCases.items():
  for k,l in j.RnaSeqFileIds.items():
    print( i , len( l.Genes ) )
# ======================================================================================================
