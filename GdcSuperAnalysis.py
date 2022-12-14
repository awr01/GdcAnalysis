from GdcLib import *
from GdcUtils import *
from scipy.stats import ttest_ind
from numpy import mean , var , log10 , warnings
from operator import itemgetter
import os , hashlib , bz2 , _pickle
warnings.filterwarnings('ignore')

# ======================================================================================================
def SortByDisease( aCase ):    
  lDiseaseType = str( aCase.DiseaseType )
  if not lDiseaseType in lCases2: lCases2[ lDiseaseType ] = []
  lCases2[ lDiseaseType ].append( aCase )
# ======================================================================================================

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
    GeneType = StarCounts.GeneCatalogue[ GeneName ].type
    if GeneType != "protein_coding" : continue
    
    mean0 , mean1  = mean( GeneData[0] ) , mean( GeneData[1] )
    if mean0 == 0 or mean1 == 0 : continue

    meanratio = mean0 / mean1    
    if ( meanratio < 3/2 ) and ( meanratio > 2/3 ) : continue   

    tscore, pvalue = ttest_ind( GeneData[0] , GeneData[1] , equal_var = False ) # Calculate the t-score between each pair of lists
    if pvalue > 1E-3 : continue
    
    logmeanratio = log10( meanratio ) / log10( 1.5 ) 
    logpvalue = -log10(pvalue)
    
    dest.write( f"{MutationOfInterest}\t{GeneName}\t{GeneType}\t{meanratio}\t{logmeanratio}\t{tscore}\t{pvalue}\t{logpvalue}\n" )
# ======================================================================================================




# ======================================================================================================
lParams = hashlib.md5( str( f"SuperAnalysis-{args.src}" ).encode('utf-8') ).hexdigest()
lFile = f".cache/{lParams}.pkl.bz2"

lCases2 = {}
# if os.path.isfile( lFile ):
  # print( "Loading cache" , flush=True )
  # with bz2.BZ2File( lFile , 'r' ) as src: lCases2 = _pickle.load( src )
# else:
LoadAndForEach( args.src , SortByDisease )
  # print( "Saving cache" , flush=True )
# with bz2.BZ2File( lFile , 'w' ) as dest: _pickle.dump( lCases2 , dest )


for DiseaseType , lCases in tqdm.tqdm( sorted( lCases2.items() ) , ncols=Ncol , desc="Analysing" ):
  lFile =  f"SuperAnalysis.{DiseaseType}.tsv"
  
  if os.path.isfile( lFile ): continue  
  
  with open( lFile , "w" ) as dest:
    dest.write( f"Mutation\tGene\tGene-type\tMut mean/WT mean\tlog_1.5(ratio)\tt-score\tp-value\tneg. log_10(p-value)\n" ) # Write headers
    
    Mutations = {}
    for lCase in lCases:
      for Gene , Mutation in lCase.Mutations.items():
        if Mutation.Classification == SilentOrSplice: continue # Ignore silents and slices
        if not Gene in Mutations: Mutations[ Gene ]  = 1
        else:                     Mutations[ Gene ] += 1  
    
    MinNumCases = 0.05 * len( lCases )
    Mutations = [ MutationOfInterest for MutationOfInterest , count in sorted( Mutations.items() , key = itemgetter(1) , reverse=True) if count >= MinNumCases ]
    
    for MutationOfInterest in tqdm.tqdm( Mutations , ncols=Ncol , leave=False , desc=DiseaseType ):
      Genes = [ ( [],[] ) for i in StarCounts.GeneCatalogue ]      
      for lCase in lCases : Analyze( lCase )  
      Finally()
# ======================================================================================================
      