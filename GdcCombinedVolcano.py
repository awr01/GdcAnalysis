from GdcLib import *
from GdcUtils import *
from openpyxl import Workbook
import matplotlib.pyplot as plt
import numpy as np

# Selected = [
#   "Brain Lower Grade Glioma",
#   "Sarcoma" ,
#   "Osteosarcoma" ,
#   "Glioblastoma Multiforme" ,
#   "Pheochromocytoma and Paraganglioma",
#   "Neuroblastoma" 
# ] 

# ======================================================================================================
def SortByDisease( aCase ):    
  index = 1                                                                      # Indices chosen for consistency with old code
  if MutationOfInterest in aCase.Mutations: 
    if aCase.Mutations[ MutationOfInterest ].Classification == SilentOrSplice : return
    index = 0

  lDiseaseType = str( aCase.DiseaseType )
  lDiseaseTypeLower = lDiseaseType.lower()

  if   "carcinoma" in lDiseaseTypeLower : key = "Carcinoma"
  elif "leukemia" in lDiseaseTypeLower:   key = "Leukemia"
  elif "melanoma" in lDiseaseTypeLower:   key = "Melanoma"
  else:                                   key = "Other"

  lCases2[ key ][ index ].append( aCase ) 
# ======================================================================================================

# ======================================================================================================
def Flatten( Data , index ):
  lRet = []      
  for j in Data:
    for i in j.StarCounts: 
      lRet.append( i.TpmUnstranded[ index ] )
  return lRet
# ======================================================================================================

# ======================================================================================================
def Finally():
  Data = {}

  # Fill the plots 
  for lDiseaseType , lCases in tqdm.tqdm( lCases2.items() , ncols=Ncol , desc="Analysing" ):

    lMut , lWt = lCases[0] , lCases[1]
    if len( lMut ) == 0 or len( lWt ) == 0 : continue
    
    x0 , y0 , x1 , y1 , x2 , y2 = [] , [] , [] , [] , [] , []
    for GeneName , Gene in tqdm.tqdm( sorted( StarCounts.GeneCatalogue.items() ) , leave=False , ncols=Ncol , desc=lDiseaseType ):
      if Gene.type != "protein_coding" : continue
      lRet = GdcStatistics( Flatten( lMut , Gene.index ) , Flatten( lWt , Gene.index ) )
      if lRet is None : continue
      if np.isnan( lRet.neg_log_pvalue ) : continue
      if ( GeneName == "DRG2" ):  
        x0.append( lRet.log_mean_ratio_with_error[0] )
        y0.append( lRet.neg_log_pvalue )
      elif ( lRet.neg_log_pvalue < 2.3 ) or ( np.fabs( lRet.log_mean_ratio_with_error[0] ) < 1 ):  
        x1.append( lRet.log_mean_ratio_with_error[0] )
        y1.append( lRet.neg_log_pvalue )
      else:
        x2.append( lRet.log_mean_ratio_with_error[0] )
        y2.append( lRet.neg_log_pvalue )   

    Data[ lDiseaseType ] = ( ( x0 , y0 ), ( x1 , y1 ), ( x2 , y2 ) )

  with bz2.BZ2File( lFile , 'w' ) as dest: _pickle.dump( Data , dest )
# ======================================================================================================

# ======================================================================================================
def DrawVolcanos():   
  fig , axs = plt.subplots( 1 , 4 , sharey=True )
  # for x in axs.flat[len(DiseaseTypeLut):]: x.set_visible( False )  
  for Cases , ax1 in zip( Data.items() , fig.axes ):
    lDiseaseType , lData = Cases
    ax1.set_xlabel( lDiseaseType , style='italic' , labelpad=1 )
    ax1.set_xlim( -25 , 25 )
    ax1.set_ylim( 1 , 200 )
    ax1.scatter( lData[1][0] , lData[1][1] , color="0.75" , s=1 )
    ax1.scatter( lData[2][0] , lData[2][1] , color="b" , s=1 )
    ax1.scatter( lData[0][0] , lData[0][1] , color="r" , s=1 )
    # ax1.text( 0.6 , 50 , f'$N_{{Mutant}} = {len(lData[0])}$\n$N_{{Wild-type}} = {len(lData[1])}$\n$p_{{value}}={lStats.pvalue:.2e}$' , fontsize="x-small" )
  plt.yscale( "log" )

  #Add the common y-axis label
  fig.add_subplot(111, frameon=False)
  plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
  plt.ylabel( "$-log_{10}($ p-value $)$" , style='italic' )
  plt.xlabel( "$log_{1.5}($ fold-ratio $)$" , style='italic', labelpad=30 )
    
  # Draw the images
  fig.set_size_inches( 16 , 20 )
  plt.tight_layout()
  fig.subplots_adjust( hspace = 0.2 , wspace = 0.0 )
  plt.savefig( 'DRG2-combined-volcano.pdf' )    
# ======================================================================================================


# ======================================================================================================
if len( args.mutations ) != 1 : raise Exception( "Exactly 1 mutation must be specified on the commandline" )
MutationOfInterest , lCases2 = args.mutations[0] , { "Carcinoma":[[],[]] , "Leukemia":[[],[]] , "Melanoma":[[],[]] , "Other":[[],[]] }         
LoadAndForEach( args.src , SortByDisease , After = Finally ) # Load src and analyze on-the-fly, making use of the optional After function
# ======================================================================================================
  