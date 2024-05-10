from GdcLib import *
from GdcUtils import *
from openpyxl import Workbook
import matplotlib.pyplot as plt
import numpy as np

# ======================================================================================================
def Classifier( aCase ):    
  DiseaseType = str( aCase.DiseaseType )
  lDiseaseType = DiseaseType.lower()

  if   "unspecified" == lDiseaseType : Class = "Unspecified"
  elif "none"        == lDiseaseType : Class = "None"
  elif "carcinoma" in lDiseaseType : Class = "Carcinoma"
  elif "leukemia" in lDiseaseType:   Class = "Leukemia"
  elif "melanoma" in lDiseaseType:   Class = "Melanoma"
  else:                              Class = "Other"

  if not Class in Dict: Dict[ Class ] = []
  if not DiseaseType in Dict[ Class ]: Dict[ Class ].append( DiseaseType )

  return Class
# ======================================================================================================


# # ======================================================================================================
# def SortByDisease( aCase ):    
#   index = 1                                                                      # Indices chosen for consistency with old code
#   if MutationOfInterest in aCase.Mutations: 
#     if aCase.Mutations[ MutationOfInterest ].Classification == SilentOrSplice : return
#     index = 0

#   lDiseaseType = str( aCase.DiseaseType )
#   if not lDiseaseType in lCases2: lCases2[ lDiseaseType ] = [ [] , [] ]
#   lCases2[ lDiseaseType ][ index ].append( aCase.CaseId )
# # ======================================================================================================

# ======================================================================================================
def Flatten( Data , index ):
  lRet = []      
  for j in Data:
    for i in j.StarCounts: 
      lRet.append( i.TpmUnstranded[ index ] )
  return lRet
# ======================================================================================================

# # ======================================================================================================
# def DataToScatter( lArgs ):
#   lDiseaseType , lCases = lArgs
#   index = multiprocessing.current_process()._identity[0]

#   if len( lCases[0] ) == 0 or len( lCases[1] ) == 0 : return

#   lMut , lWt = LoadCases( args.src , lCases[0] ) , LoadCases( args.src , lCases[1] )
  
#   x0 , y0 , x1 , y1 , x2 , y2 , x3 , y3 , x4 , y4 = [] , [] , [] , [] , [] , [] , [] , [] , [] , []
#   for GeneName , Gene in tqdm.tqdm( sorted( StarCounts.GeneCatalogue.items() ) , leave=False , ncols=Ncol , desc=lDiseaseType, position=index ):
#     if Gene.type != "protein_coding" : continue
#     lRet = GdcStatistics( Flatten( lMut , Gene.index ) , Flatten( lWt , Gene.index ) )
#     if lRet is None : continue
#     if np.isnan( lRet.neg_log_pvalue ) : continue
#     if ( GeneName == "DRG2" ):  
#       x0.append( lRet.log_mean_ratio_with_error[0] )
#       y0.append( lRet.neg_log_pvalue )
#     elif ( GeneName == "TRIM24" ):  
#       x3.append( lRet.log_mean_ratio_with_error[0] )
#       y3.append( lRet.neg_log_pvalue )
#     elif ( GeneName == "SLFN11" ):  
#       x4.append( lRet.log_mean_ratio_with_error[0] )
#       y4.append( lRet.neg_log_pvalue )      
#     elif ( lRet.neg_log_pvalue < 2.3 ) or ( np.fabs( lRet.log_mean_ratio_with_error[0] ) < 1 ):  
#       x1.append( lRet.log_mean_ratio_with_error[0] )
#       y1.append( lRet.neg_log_pvalue )
#     else:
#       x2.append( lRet.log_mean_ratio_with_error[0] )
#       y2.append( lRet.neg_log_pvalue )   

#   return lDiseaseType , ( ( x0 , y0 ), ( x1 , y1 ), ( x2 , y2 ), ( x3 , y3 ), ( x4 , y4 ) )  
# # ======================================================================================================


# ======================================================================================================
def ForEachClass( Class , Cases , index ):
  lMut , lWt = [] , []

  for Case in Cases:
    if MutationOfInterest in Case.Mutations: 
      if Case.Mutations[ MutationOfInterest ].Classification != SilentOrSplice : lMut.append( Case )   
    else: lWt.append( Case )  

  if len( lMut ) == 0 or len( lWt ) == 0 : return
  
  x0 , y0 , x1 , y1 , x2 , y2 = [] , [] , [] , [] , [] , []
  for GeneName , Gene in tqdm.tqdm( sorted( StarCounts.GeneCatalogue.items() ) , leave=False , ncols=Ncol , desc=f"{Class}: Analysing", position=index ):
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

  return ( ( x0 , y0 ), ( x1 , y1 ), ( x2 , y2 ) )  
# ======================================================================================================


# # ======================================================================================================
# def Finally():
#   # Fill the plots 
#   with multiprocessing.Pool() as pool:
#     for res in tqdm.tqdm( pool.imap( DataToScatter , sorted( lCases2.items() ) ) , ncols=Ncol , desc="Analysing" ):
#       if not res is None: Data[ res[0] ] = res[1]
# # ======================================================================================================



# ======================================================================================================
def DrawVolcanos( Data ):    
  # Draw the plots
  fig , axs = plt.subplots( 2 , 3 , sharey=True )
  for x in axs.flat[len(Data):]: x.set_visible( False )  
  for (lDiseaseType , lData) , ax1 in zip( Data.items() , fig.axes ):
    if lData is None: continue

    ax1.set_xlabel( lDiseaseType , style='italic' , labelpad=1 )
    ax1.set_xlim( -25 , 25 )
    ax1.set_ylim( 1/200 , 200 )
    ax1.scatter( lData[1][0] , lData[1][1] , color="0.75" , s=1 )
    ax1.scatter( lData[2][0] , lData[2][1] , color="b" , s=1 )
    ax1.scatter( lData[0][0] , lData[0][1] , color="r" , s=1 )

    ax1.text( -24 , 0.01 , "\n".join( sorted( Dict[lDiseaseType] ) ) , fontsize="x-small" )
    ax1.grid( True )

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
MutationOfInterest, Dict = args.mutations[0], {}

# lFile = f".cache/GdcVolcanoRaw.pkl.bz"

# if os.path.isfile( lFile ):
#   print( "> Using cached stats" , flush=True )
#   with bz2.BZ2File( lFile , 'r' ) as src: Data = _pickle.load( src )  
# else:
#   print( "> Extracting stats" , flush=True )
# LoadAndForEach( args.src , SortByDisease , After = Finally ) # Load src and analyze on-the-fly, making use of the optional After function
# with bz2.BZ2File( lFile , 'w' ) as dest: _pickle.dump( Data , dest )

LoadAndClassify( args.src , Classifier , ForEachClass , DrawVolcanos , maxthreads=3 )
# ======================================================================================================
  