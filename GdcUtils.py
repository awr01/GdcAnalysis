import argparse
from GdcUtils import *

# Genes = [ "ATRX" , "GAPDH" , "TUBA1A" , "ACTB" , "ALB" , "ALOX12" , "ANGPTL7" , "AOX1" , "APOE" , "ATOX1" , 
          # "BNIP3" , "CAT" , "CCL5" , "CCS" , "CSDE1" , "CYBA" , "CYGB" , "DGKK" , "DHCR24" , "DUOX1" , "DUOX2" , 
          # "DUSP1" , "EPHX2" , "EPX" , "FOXM1" , "GLRX2" , "GPR156" , "GPX1" , "GPX2" , "GPX3" , "GPX4" , "GPX5" , 
          # "GPX6" , "GPX7" , "GSR" , "GSS" , "GSTZ1" , "GTF2I" , "KRT1" , "LPO" , "MBL2" , "MGST3" , "MPO" , 
          # "MPV17" , "MSRA" , "MT3" , "TESMIN" , "NCF1" , "NCF1B" , "NCF1C" , "NCF2" , "NME5" , "NOS2" , "NOX5" , 
          # "NUDT1" , "OXR1" , "OXSR1" , "PDLIM1" , "IPCEF1" , "PNKP" , "PRDX1" , "PRDX2" , "PRDX3" , "PRDX4" , "PRDX5" , "PRDX6" , 
          # "PREX1" , "PRG3" , "PRNP" , "PTGS1" , "PTGS2" , "PXDN" , "PXDNL" , "RNF7" , "SCARA3" , "SELENOS" , "SELENOP" , "SFTPD" , 
          # "SGK2" , "SIRT2" , "SOD1" , "SOD2" , "SOD3" , "SRXN1" , "STK25" , "TPO" , "TTN" , "TXNDC2" , "TXNRD1" , "TXNRD2" ]

# Mutations = [ "ATRX" , "IDH1" , "SETD2" ]

parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')

parser.add_argument( '--genes' , nargs='+' , default=[] , help='Genes of interest')
parser.add_argument( '--mutations' , nargs='+' , default=[] , help='Mutations of interest')

parser.add_argument( '--gene-file' , help='Genes of interest')
parser.add_argument( '--mutation-file' , help='Mutations of interest')

args = parser.parse_args()

if not args.src.endswith( ".tgz" ): raise Exception( "Source file must have '.tgz' file-extension" )
if not args.dest.endswith( ".tsv" ): raise Exception( "Destination file must have '.tsv' file-extension" )

if not args.gene_file is None:
  with open( args.gene_file , "r" ) as src:
    for line in src: args.genes.append( line.strip() )

if not args.mutation_file is None:
  with open( args.mutation_file , "r" ) as src:
    for line in src: args.mutations.append( line.strip() )
