import argparse
from GdcLib import *

# ======================================================================================================
def Init():
  dest.write( "CaseId" )                                                         # Write the header-line to the destination tsv
  for i in Mutations : dest.write( "\tMutations-" + i )
  for i in Genes :     dest.write( "\tStarCount-" + i )
  dest.write( "\n" )
# ======================================================================================================

# ======================================================================================================
def Analyze( aCase ):
  for lStarCount in aCase.StarCounts:                                            # Then iterate over each star-count file
    dest.write( aCase.CaseId )
    for i in Mutations: dest.write( "\t" + str( aCase.GetMutations( i , "" ) ) )         
    for i in Genes:     dest.write( "\t" + str( lStarCount[i] ) )
    dest.write( "\n" ) 
# ======================================================================================================



# ======================================================================================================
Genes = [ "ATRX" , "GAPDH" , "TUBA1A" , "ACTB" , "ALB" , "ALOX12" , "ANGPTL7" , "AOX1" , "APOE" , "ATOX1" , 
          "BNIP3" , "CAT" , "CCL5" , "CCS" , "CSDE1" , "CYBA" , "CYGB" , "DGKK" , "DHCR24" , "DUOX1" , "DUOX2" , 
          "DUSP1" , "EPHX2" , "EPX" , "FOXM1" , "GLRX2" , "GPR156" , "GPX1" , "GPX2" , "GPX3" , "GPX4" , "GPX5" , 
          "GPX6" , "GPX7" , "GSR" , "GSS" , "GSTZ1" , "GTF2I" , "KRT1" , "LPO" , "MBL2" , "MGST3" , "MPO" , 
          "MPV17" , "MSRA" , "MT3" , "TESMIN" , "NCF1" , "NCF1B" , "NCF1C" , "NCF2" , "NME5" , "NOS2" , "NOX5" , 
          "NUDT1" , "OXR1" , "OXSR1" , "PDLIM1" , "IPCEF1" , "PNKP" , "PRDX1" , "PRDX2" , "PRDX3" , "PRDX4" , "PRDX5" , "PRDX6" , 
          "PREX1" , "PRG3" , "PRNP" , "PTGS1" , "PTGS2" , "PXDN" , "PXDNL" , "RNF7" , "SCARA3" , "SELENOS" , "SELENOP" , "SFTPD" , 
          "SGK2" , "SIRT2" , "SOD1" , "SOD2" , "SOD3" , "SRXN1" , "STK25" , "TPO" , "TTN" , "TXNDC2" , "TXNRD1" , "TXNRD2" ]

Mutations = [ "ATRX" , "IDH1" , "SETD2" ]

parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')
args = parser.parse_args()
if not args.src.endswith( ".tgz" ): raise Exception( "Source file must have '.tgz' file-extension" )
if not args.dest.endswith( ".tsv" ): raise Exception( "Destination file must have '.tsv' file-extension" )

with open( args.dest , "w" ) as dest: LoadAndForEach( args.src , Analyze , Init ) # Open the destination tsv and call the file-handle 'dest', then load src and analyze on-the-fly, making use of the optional Before argument
# ======================================================================================================
