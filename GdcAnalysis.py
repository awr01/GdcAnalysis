from GdcLib import *
from GdcUtils import *
from scipy.stats import ttest_ind
from numpy import mean , var , log10 , warnings
warnings.filterwarnings('ignore')

# ======================================================================================================
def Init():
  global Genes
  Genes = [ ( [],[] ) for i in StarCounts.GeneCatalogue ]
  dest.write( f"Gene\tFlagged\tMut-count\tMut-mean\tMut-var\tWT-count\tWT-mean\tWT-var\tMut mean/WT mean\tlog_1.5(ratio)\tt-score\tp-value\t'-log_10(p-value)\n" ) # Write headers
# ======================================================================================================

# ======================================================================================================
def Analyze( aCase ):    
  index = 1                                                                      # Indices chosen for consistency with old code
  if MutationOfInterest in aCase.Mutations: 
    if aCase.Mutations[ MutationOfInterest ].Classification == SilentOrSplice : return       # Ignore silents and slices
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
    logmeanratio = log10( meanratio ) / log10( 1.5 ) 
    
    flag = "*" if k in args.genes else ""
    t01, p01 = ttest_ind( v[0] , v[1] , equal_var = False ) # Calculate the t-score between each pair of lists

    dest.write( f"{k}\t{flag}\t{num0}\t{mean0}\t{var0}\t{num1}\t{mean1}\t{var1}\t{meanratio}\t{logmeanratio}\t{t01}\t{p01}\t{-log10(p01)}\n" )
# ======================================================================================================

# ======================================================================================================
if len( args.mutations ) != 1 : raise Exception( "Exactly 1 mutation must be specified on the commandline" )
MutationOfInterest , Genes = args.mutations[0] , None         
with open( args.dest , "w" ) as dest: LoadAndForEach( args.src , Analyze , Init , Finally ) # Open the destination tsv and call the file-handle 'dest', then load src and analyze on-the-fly, making use of the optional Before and After functions 
# ======================================================================================================
