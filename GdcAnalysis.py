import argparse
from GdcLib import *
from scipy.stats import ttest_ind
from numpy import mean , var , log10 , warnings
warnings.filterwarnings('ignore')

log_scale = 1.0 / log10( 1.5 ) 

# ======================================================================================================
def Init():
  global Genes
  if Genes is None: Genes = [ ( [],[] ) for i in StarCounts.GeneCatalogue ]
  dest.write( f"Gene\tFlagged\tMut-count\tMut-mean\tMut-var\tWT-count\tWT-mean\tWT-var\tMut mean/WT mean\tlog_1.5(ratio)\tt-score\tp-value\t'-log_10(p-value)\n" ) # Write headers
# ======================================================================================================

# ======================================================================================================
def Analyze( aCase ):    
  index = 1                                                                      # Indices chosen for consistency with old code
  if "ATRX" in aCase.Mutations: 
    if aCase.Mutations[ "ATRX" ].Classification == SilentOrSplice : return       # Ignore silents and slices
    index = 0
        
  for lStarCount in aCase.StarCounts:                                            # Then iterate over each star-count file
    for i , j in zip( lStarCount.Genes , Genes ):
      if not i is None: j[ index ].append( i ) 
# ======================================================================================================

# ======================================================================================================
def Finally():
  for k , v in tqdm.tqdm( zip( StarCounts.GeneCatalogue , Genes ) , ncols=Ncol , total=len(Genes) , desc="Saving" ):
    num0 , mean0 , var0 = len( v[0] ) , mean( v[0] ) , var( v[0] )               # Calculate the mean and variance of the muts
    num1 , mean1 , var1 = len( v[1] ) , mean( v[1] ) , var( v[1] )               # Calculate the mean and variance of the WTs
    if mean0 == 0 or mean1 == 0 : continue

    meanratio = mean0 / mean1
    logmeanratio = log10( meanratio ) * log_scale
    
    flag = "*" if k in GenesOfInterest else ""
    t01, p01 = ttest_ind( v[0] , v[1] , equal_var = False ) # Calculate the t-score between each pair of lists

    dest.write( f"{k}\t{flag}\t{num0}\t{mean0}\t{var0}\t{num1}\t{mean1}\t{var1}\t{meanratio}\t{logmeanratio}\t{t01}\t{p01}\t{-log10(p01)}\n" )
# ======================================================================================================



# ======================================================================================================
GenesOfInterest = [ "ATRX" , "GAPDH" , "TUBA1A" , "ACTB" , "ALB" , "ALOX12" , "ANGPTL7" , "AOX1" , "APOE" , "ATOX1" , 
          "BNIP3" , "CAT" , "CCL5" , "CCS" , "CSDE1" , "CYBA" , "CYGB" , "DGKK" , "DHCR24" , "DUOX1" , "DUOX2" , 
          "DUSP1" , "EPHX2" , "EPX" , "FOXM1" , "GLRX2" , "GPR156" , "GPX1" , "GPX2" , "GPX3" , "GPX4" , "GPX5" , 
          "GPX6" , "GPX7" , "GSR" , "GSS" , "GSTZ1" , "GTF2I" , "KRT1" , "LPO" , "MBL2" , "MGST3" , "MPO" , 
          "MPV17" , "MSRA" , "MT3" , "TESMIN" , "NCF1" , "NCF1B" , "NCF1C" , "NCF2" , "NME5" , "NOS2" , "NOX5" , 
          "NUDT1" , "OXR1" , "OXSR1" , "PDLIM1" , "IPCEF1" , "PNKP" , "PRDX1" , "PRDX2" , "PRDX3" , "PRDX4" , "PRDX5" , "PRDX6" , 
          "PREX1" , "PRG3" , "PRNP" , "PTGS1" , "PTGS2" , "PXDN" , "PXDNL" , "RNF7" , "SCARA3" , "SELENOS" , "SELENOP" , "SFTPD" , 
          "SGK2" , "SIRT2" , "SOD1" , "SOD2" , "SOD3" , "SRXN1" , "STK25" , "TPO" , "TTN" , "TXNDC2" , "TXNRD1" , "TXNRD2" ]

parser = argparse.ArgumentParser() # Read the commandline arguments
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')
args = parser.parse_args()
if not args.src.endswith( ".tgz" ): raise Exception( "Source file must have '.tgz' file-extension" )
if not args.dest.endswith( ".tsv" ): raise Exception( "Destination file must have '.tsv' file-extension" )
          
with open( args.dest , "w" ) as dest:                                       # Open the destination tsv and call the file-handle 'dest'
  Genes = None
  LoadAndForEach( args.src , Analyze , Init , Finally )  # Load the formatted tarball and analyze on-the-fly, making use of the optional Before and After functions 
# ======================================================================================================
