from GdcLib import *
from GdcUtils import *
from scipy.stats import ttest_ind
from numpy import mean , var , log10 , warnings
from operator import itemgetter
warnings.filterwarnings('ignore')

# ======================================================================================================
def Analyze( aCase ):    
  index = 1                                                                      # Indices chosen for consistency with old code
  if MutationOfInterest in aCase.Mutations: 
    if aCase.Mutations[ MutationOfInterest ].Classification == SilentOrSplice : return       # Ignore silents and slices
    index = 0
        
  for lStarCount in aCase.StarCounts:                                            # Then iterate over each star-count file
    for TpmUnstranded , GeneData in zip( lStarCount.TpmUnstranded , Genes ):
      if not TpmUnstranded is None: GeneData[ index ].append( TpmUnstranded ) 
# ======================================================================================================

# ======================================================================================================
def Finally():
  for GeneName , GeneData in zip( StarCounts.GeneCatalogue , Genes ):  
    mean0 , mean1  = mean( GeneData[0] ) , mean( GeneData[1] )
    if mean0 == 0 or mean1 == 0 : continue

    meanratio = mean0 / mean1    
    if ( meanratio < 3/2 ) and ( meanratio > 2/3 ) : continue   

    t01, p01 = ttest_ind( GeneData[0] , GeneData[1] , equal_var = False ) # Calculate the t-score between each pair of lists
    if p01 > 1E-30 : continue

    Results.append( ( MutationOfInterest , GeneName , meanratio , t01 , p01 ) )
# ======================================================================================================

# ======================================================================================================
def Save():
  dest.write( f"Mutation\tGene\tGene-type\tMut mean/WT mean\tlog_1.5(ratio)\tt-score\tp-value\tneg. log_10(p-value)\n" ) # Write headers

  for Result in sorted( Results , key = itemgetter(4) ):
    Mutation , GeneName , meanratio , t01 , p01 = Result    
    logmeanratio = log10( meanratio ) / log10( 1.5 ) 
    GeneType = StarCounts.GeneCatalogue[ GeneName ].type
    dest.write( f"{Mutation}\t{GeneName}\t{GeneType}\t{meanratio}\t{logmeanratio}\t{t01}\t{p01}\t{-log10(p01)}\n" )
# ======================================================================================================

# ======================================================================================================
Genes = None
Results = []
         
with open( args.dest , "w" ) as dest:
  lCases = LoadCases( args.src )
  
  Mutations = {}
  for lCase in lCases:
    for Gene , Mutation in lCase.Mutations.items():
      if Mutation.Classification == SilentOrSplice: continue # Ignore silents and slices
      if not Gene in Mutations: Mutations[ Gene ]  = 1
      else:                     Mutations[ Gene ] += 1

  MinNumCases = 0.05 * len( lCases )
  Mutations = [ MutationOfInterest for MutationOfInterest , count in sorted( Mutations.items() , key = itemgetter(1) , reverse=True) if count >= MinNumCases ]

  for MutationOfInterest in tqdm.tqdm( Mutations , ncols=Ncol , desc="Scanning mutations" ):
    Genes = [ ( [],[] ) for i in StarCounts.GeneCatalogue ]      
    for lCase in lCases : Analyze( lCase )  
    Finally()

  Save()
# ======================================================================================================
