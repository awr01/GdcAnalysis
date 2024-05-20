# import argparse
from GdcLib import *
# import os , hashlib , bz2 , _pickle , math , gzip


# parser = argparse.ArgumentParser()
# parser.add_argument( '--src' , required=True , help='The source tarball')
# parser.add_argument( '--dest' , required=True , help='The destination file')

# parser.add_argument( '--mutation' , required=True , help='Mutation of interest')

# # parser.add_argument( '--genes' , nargs='+' , default=[] , help='Genes of interest')
# # parser.add_argument( '--mutations' , nargs='+' , default=[] , help='Mutations of interest')

# # parser.add_argument( '--gene-file' , help='Genes of interest')
# # parser.add_argument( '--mutation-file' , help='Mutations of interest')

# args = parser.parse_args()

# if not args.src .endswith( ".tar"  ): raise Exception( "Source file must have '.tar' file-extension" )
# # if not args.dest.endswith( ".xlsx" ): raise Exception( "Destination file must have '.xlsx' file-extension" )

 

# # if not args.gene_file is None:
# #   with open( args.gene_file , "r" ) as src:
# #     for line in src: args.genes.append( line.strip() )

# # if not args.mutation_file is None:
# #   with open( args.mutation_file , "r" ) as src:
# #     for line in src: args.mutations.append( line.strip() )





# ======================================================================================================
from scipy.stats import ttest_ind , sem
from numpy import mean , std , sqrt , log10 , log , warnings
warnings.filterwarnings( 'ignore' )

const0 = 1.0 / log10( 1.5 )
const1 = const0 / log( 10.0 )
# ======================================================================================================


# ======================================================================================================
class GdcStatsStruct:
  class GdcStatsStructInner:
    def __init__( self ):
      self.count = None
      self.mean = None
      self.mean_error = None
      self.sd = None

  def __init__( self ):
    self.mut = GdcStatsStruct.GdcStatsStructInner()
    self.wt  = GdcStatsStruct.GdcStatsStructInner()

    self.mean_ratio_with_error     = None
    self.log_mean_ratio_with_error = None

    self.tscore = None
    self.pvalue = None
    self.neg_log_pvalue = None
# ======================================================================================================


# ======================================================================================================
def GdcStatistics( aMut , aWT , aRatioCut = False , aPvalueCut = None ):  
  lRet = GdcStatsStruct()
  
  lRet.mut.count  = len( aMut )
  if lRet.mut.count == 0 : return None
  
  lRet.wt.count = len( aWT )
  if lRet.wt.count == 0 : return None

  lRet.mut.mean = mean( aMut ) 
  if lRet.mut.mean == 0 : return None
  
  lRet.wt.mean = mean( aWT )
  if lRet.wt.mean == 0 : return None

  mean_ratio = lRet.mut.mean / lRet.wt.mean  

  if aRatioCut:
    if ( mean_ratio < 3/2 ) and ( mean_ratio > 2/3 ) : return None 
    
  lRet.tscore , lRet.pvalue = ttest_ind( aMut , aWT , equal_var = False ) # Calculate the t-score between each pair of lists

  if not aPvalueCut is None:
    if lRet.pvalue > aPvalueCut : return None

  lRet.mut.mean_error , lRet.mut.sd = sem( aMut ) , std( aMut )
  lRet.wt.mean_error  , lRet.wt.sd  = sem( aWT )  , std( aWT )

  mut_rel_error , wt_rel_error      = lRet.mut.mean_error/lRet.mut.mean , lRet.wt.mean_error/lRet.wt.mean
  mean_ratio_rel_error              = sqrt( (mut_rel_error*mut_rel_error) + (wt_rel_error*wt_rel_error) )
  
  lRet.mean_ratio_with_error        = ( mean_ratio   , mean_ratio_rel_error * mean_ratio )
  lRet.log_mean_ratio_with_error    = ( log10( mean_ratio ) * const0 , mean_ratio_rel_error * const1 )  
  lRet.neg_log_pvalue               = -log10( lRet.pvalue )  

  return lRet
# ======================================================================================================


# ======================================================================================================
def FlattenTpmUnstranded( Data , index ):
  lRet = []      
  for j in Data:
    for i in j.StarCounts: 
      lRet.append( i.TpmUnstranded[ index ] )
  return lRet
# ======================================================================================================

# ======================================================================================================
def SeparateMutandAndWildType( Cases , Mutation ):
  lMut , lWt , Diseases = [] , [] , {}

  for Case in Cases:

    info = str( Case.DiseaseType )
    if not info in Diseases: Diseases[ info ] = { "Mut":0 , "WT":0 }

    if Mutation in Case.Mutations: 
      if Case.Mutations[ Mutation ].Classification == SilentOrSplice : continue
      lMut.append( Case ) 
      Diseases[ info ][ "Mut" ] += 1  
    else:
      lWt.append( Case )  
      Diseases[ info ][ "WT" ] += 1  

  return lMut , lWt , Diseases
# ======================================================================================================



# # ======================================================================================================
# def GdcAnalysis_ForEachClass( Class , Cases , index ):
#   lMut , lWt , Diseases = SeparateMutandAndWildType( Cases )
#   if len( lMut ) == 0 or len( lWt ) == 0 : return
  
#   Results = {}
#   for GeneName , Gene in tqdm.tqdm( sorted( StarCounts.GeneCatalogue.items() ) , leave=False , ncols=Ncol , desc=f"{Class}: Analysing", position=index ):
#     if Gene.type != "protein_coding" : continue
#     lRet = GdcStatistics( FlattenTpmUnstranded( lMut , Gene.index ) , FlattenTpmUnstranded( lWt , Gene.index ) )
#     if lRet is None: continue
#     if math.isnan( lRet.neg_log_pvalue ) : continue
#     Results[ GeneName ] = lRet
#   return Diseases , Results

# def GdcAnalysis( ClassifyFn , ClassifyDesc , FinallyFn , maxthreads=None ): 
#   CacheFile = f".cache/GdcAnalysis.{args.mutation}.{ClassifyDesc}.pkl.gz"
#   LoadAndClassify( args.src , ClassifyFn , GdcAnalysis_ForEachClass , FinallyFn , maxthreads=maxthreads , cachefile=CacheFile )    
# # ======================================================================================================
