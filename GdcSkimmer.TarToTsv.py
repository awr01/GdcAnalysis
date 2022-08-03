import argparse
from GdcLib import *

parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')
args = parser.parse_args()

# ==============================================================================================================================================
# Extract the list of genes to a single tsv    

Genes = [ "ATRX" , "GAPDH" , "TUBA1A" , "ACTB" , "ALB" , "ALOX12" , "ANGPTL7" , "AOX1" , "APOE" , "ATOX1" , 
          "BNIP3" , "CAT" , "CCL5" , "CCS" , "CSDE1" , "CYBA" , "CYGB" , "DGKK" , "DHCR24" , "DUOX1" , "DUOX2" , 
          "DUSP1" , "EPHX2" , "EPX" , "FOXM1" , "GLRX2" , "GPR156" , "GPX1" , "GPX2" , "GPX3" , "GPX4" , "GPX5" , 
          "GPX6" , "GPX7" , "GSR" , "GSS" , "GSTZ1" , "GTF2I" , "KRT1" , "LPO" , "MBL2" , "MGST3" , "MPO" , 
          "MPV17" , "MSRA" , "MT3" , "TESMIN" , "NCF1" , "NCF1B" , "NCF1C" , "NCF2" , "NME5" , "NOS2" , "NOX5" , 
          "NUDT1" , "OXR1" , "OXSR1" , "PDLIM1" , "IPCEF1" , "PNKP" , "PRDX1" , "PRDX2" , "PRDX3" , "PRDX4" , "PRDX5" , "PRDX6" , 
          "PREX1" , "PRG3" , "PRNP" , "PTGS1" , "PTGS2" , "PXDN" , "PXDNL" , "RNF7" , "SCARA3" , "SELENOS" , "SELENOP" , "SFTPD" , 
          "SGK2" , "SIRT2" , "SOD1" , "SOD2" , "SOD3" , "SRXN1" , "STK25" , "TPO" , "TTN" , "TXNDC2" , "TXNRD1" , "TXNRD2" ]

with open( args.dest , "w" ) as dest:                                    # Open the destination tsv and call the file-handle 'dest'
  lCases = LoadFormattedTarball( args.src , Genes )                      # Load the GDC data from the formatted tarball as a dictionary of CaseId strings and Case objects
                                                                         # keeping only the specified Genes from the star-count file

  dest.write( "CaseId\tATRX-mutation\t" + "\t".join( Genes ) + "\n" )    # Write the header-line to the destination tsv 
  for lCaseId , lCase in sorted( lCases.items() ):                       # Sort the cases by case-id, then iterate over them storing the CaseId and Case in separate variables

    lMutations = lCase.getMutations( "ATRX" )                            # Get any mutations of ATRX
    if lMutations is None : lMutations = "No WXS Data"                   # If there is no info, say so
    else: lMutations = " , ".join( [ str(i) for i in lMutations ] )      # Else pack all mutations into a single comma-delimited string
    
    if lCase.StarCounts is None:                                         # If there is no star-count info...
      dest.write( f"{lCaseId}\t{lMutations}\tNo star-count\n" )          #   Write the information that is available to the destination file using python 3's f-strings
    else:
      for lStarCount in lCase.StarCounts:                                # Else, iterate over each star-count file
        lStarCount = [ str( lStarCount[i] ) for i in Genes ]             #   Convert the dictionary for that star-count file into a list of strings in the same order as the genes list
        lStarCount = "\t".join( lStarCount )                             #   Convert the list to a single tab-delimeted string
        dest.write( f"{lCaseId}\t{lMutations}\t{lStarCount}\n" )         #   Write the information to the destination file using python 3's f-strings
# ==============================================================================================================================================

