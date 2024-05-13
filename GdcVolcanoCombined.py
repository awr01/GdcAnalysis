from GdcLib import *
from GdcUtils import *
import matplotlib.pyplot as plt
import numpy as np

# ======================================================================================================
def Classifier( aCase ):    
  DiseaseType = str( aCase.DiseaseType )
  lDiseaseType = DiseaseType.lower()

  if "unspecified" == lDiseaseType: Class = "Unspecified"
  elif "none"      == lDiseaseType: Class = "None"
  elif "carcinoma" in lDiseaseType: Class = "Carcinoma"
  elif "leukemia"  in lDiseaseType: Class = "Leukemia"
  elif "melanoma"  in lDiseaseType: Class = "Melanoma"
  else:                             Class = "Other"

  if not Class in GdcStatsStruct.Aux.Dict: GdcStatsStruct.Aux.Dict[ Class ] = []
  if not DiseaseType in GdcStatsStruct.Aux.Dict[ Class ]: GdcStatsStruct.Aux.Dict[ Class ].append( DiseaseType )

  return Class
# ======================================================================================================


# # ======================================================================================================
# def Flatten( Data , index ):
#   lRet = []      
#   for j in Data:
#     for i in j.StarCounts: 
#       lRet.append( i.TpmUnstranded[ index ] )
#   return lRet

# def ForEachClass( Class , Cases , index ):
#   lMut , lWt = [] , []

#   for Case in Cases:
#     if MutationOfInterest in Case.Mutations: 
#       if Case.Mutations[ MutationOfInterest ].Classification != SilentOrSplice : lMut.append( Case )   
#     else: lWt.append( Case )  

#   if len( lMut ) == 0 or len( lWt ) == 0 : return
  
#   x0 , y0 , x1 , y1 , x2 , y2 = [] , [] , [] , [] , [] , []
#   for GeneName , Gene in tqdm.tqdm( sorted( StarCounts.GeneCatalogue.items() ) , leave=False , ncols=Ncol , desc=f"{Class}: Analysing", position=index ):
#     if Gene.type != "protein_coding" : continue
#     lRet = GdcStatistics( Flatten( lMut , Gene.index ) , Flatten( lWt , Gene.index ) )
#     if lRet is None : continue
#     if np.isnan( lRet.neg_log_pvalue ) : continue
#     if ( GeneName == "DRG2" ):  
#       x0.append( lRet.log_mean_ratio_with_error[0] )
#       y0.append( lRet.neg_log_pvalue )   
#     elif ( lRet.neg_log_pvalue < 2.3 ) or ( np.fabs( lRet.log_mean_ratio_with_error[0] ) < 1 ):  
#       x1.append( lRet.log_mean_ratio_with_error[0] )
#       y1.append( lRet.neg_log_pvalue )
#     else:
#       x2.append( lRet.log_mean_ratio_with_error[0] )
#       y2.append( lRet.neg_log_pvalue )   

#   return ( ( x0 , y0 ), ( x1 , y1 ), ( x2 , y2 ) )  
# # ======================================================================================================


# # ======================================================================================================
# def DrawVolcanos( Data ):    
#   # Draw the plots
#   fig , axs = plt.subplots( 2 , 3 , sharey=True )
#   for x in axs.flat[len(Data):]: x.set_visible( False )  
#   for (lDiseaseType , lData) , ax1 in zip( Data.items() , fig.axes ):
#     if lData is None: continue

#     ax1.set_xlabel( lDiseaseType , style='italic' , labelpad=1 )
#     ax1.set_xlim( -25 , 25 )
#     ax1.set_ylim( 1/200 , 200 )
#     ax1.scatter( lData[1][0] , lData[1][1] , color="0.75" , s=1 )
#     ax1.scatter( lData[2][0] , lData[2][1] , color="b" , s=1 )
#     ax1.scatter( lData[0][0] , lData[0][1] , color="r" , s=1 )
#     ax1.grid( True )
#     ax1.text( -24 , 0.01 , "\n".join( sorted( GdcStatsStruct.Aux.Dict[lDiseaseType] ) ) , fontsize="x-small" )

#   plt.yscale( "log" )

#   #Add the common y-axis label
#   fig.add_subplot(111, frameon=False)
#   plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
#   plt.ylabel( "$-log_{10}($ p-value $)$" , style='italic' )
#   plt.xlabel( "$log_{1.5}($ fold-ratio $)$" , style='italic', labelpad=30 )
    
#   # Draw the images
#   fig.set_size_inches( 16 , 20 )
#   plt.tight_layout()
#   fig.subplots_adjust( hspace = 0.2 , wspace = 0.0 )
#   plt.savefig( 'DRG2-combined-volcano.pdf' )    
# # ======================================================================================================


# # ======================================================================================================
# if len( args.mutations ) != 1 : raise Exception( "Exactly 1 mutation must be specified on the commandline" )
# MutationOfInterest, GdcStatsStruct.Aux.Dict = args.mutations[0], {}
# LoadAndClassify( args.src , Classifier , ForEachClass , DrawVolcanos , maxthreads=3 )
# # ======================================================================================================
  

# ======================================================================================================
def DrawVolcanos( Data ):    
  Data = { k:v for k,v in Data.items() if (not v is None) and len(v) }

  n = len(Data)
  a = int( np.ceil( np.sqrt( n ) ) )
  b = int( np.ceil( n/a ) )
  fig , axs = plt.subplots( b , a , sharey=True )
  for x in axs.flat[n:]: x.set_visible( False )

  for (lDiseaseType , lData) , ax1 in tqdm.tqdm( list( zip( Data.items() , fig.axes ) ) , ncols=Ncol , desc=f"Drawing Volcanos" ):

    x0 , y0 , x1 , y1 , x2 , y2 = [] , [] , [] , [] , [] , []

    for GeneName, Stats in lData.items():
      if ( GeneName == "DRG2" ):  
        x0.append( Stats.log_mean_ratio_with_error[0] )
        y0.append( Stats.neg_log_pvalue )   
      elif ( Stats.neg_log_pvalue < 2.3 ) or ( np.fabs( Stats.log_mean_ratio_with_error[0] ) < 1 ):  
        x1.append( Stats.log_mean_ratio_with_error[0] )
        y1.append( Stats.neg_log_pvalue )
      else:
        x2.append( Stats.log_mean_ratio_with_error[0] )
        y2.append( Stats.neg_log_pvalue )   

    # ax1.set_xlabel( lDiseaseType , style='italic' , labelpad=1 )
    ax1.set_xlim( -25 , 25 )
    ax1.set_ylim( 1/200 , 200 )
    ax1.scatter( x1 , y1 , color="0.75" , s=1 )
    ax1.scatter( x2 , y2 , color="b" , s=1 )
    ax1.scatter( x0 , y0 , color="r" , s=1 )
    ax1.grid( True )
    ax1.text( -24 , 0.01 , "\n".join( sorted( GdcStatsStruct.Aux.Dict[lDiseaseType] ) ) , fontsize="x-small" )

    # if len( lDiseaseType ) > 40 :
    #   index = lDiseaseType.rfind( ' ' , 0 , 40 )
    #   lDiseaseType = lDiseaseType[:index] + '\n' + lDiseaseType[index:]

    # ax1.text( -24 , 0.01 , lDiseaseType , fontsize="x-small" )

  #Add the common y-axis label
  plt.yscale( "log" )
  fig.add_subplot(111, frameon=False)
  plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
  plt.ylabel( "$-log_{10}($ p-value $)$" , style='italic' )
  plt.xlabel( "$log_{1.5}($ fold-ratio $)$" , style='italic', labelpad=30 )
    
  # Draw the images
  fig.set_size_inches( 16 , 20 )
  plt.tight_layout()
  fig.subplots_adjust( hspace = 0.2 , wspace = 0.0 )
  plt.savefig( args.dest )    
# ======================================================================================================


# ====================================================================================================== 
if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
setattr( GdcStatsStruct.Aux , "Dict" , {} )
GdcAnalysis( Classifier , "BroadClasses" , DrawVolcanos , maxthreads=3 )
# ======================================================================================================  