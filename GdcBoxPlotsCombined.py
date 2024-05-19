from GdcLib import *
from GdcUtils import *
from openpyxl import Workbook
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

  return Class
# ======================================================================================================

# ======================================================================================================
def BoxPlot_ForEachClass( Class , Cases , index ):
  lMut , lWt , Diseases = SeparateMutandAndWildType( Cases )
  if len( lMut ) == 0 or len( lWt ) == 0 : return

  DRG2 = StarCounts.GeneCatalogue[ "DRG2" ].index    
  return Diseases, ( FlattenTpmUnstranded( lMut , DRG2 ) , FlattenTpmUnstranded( lWt , DRG2 ) )
# ======================================================================================================

# ======================================================================================================
def DrawBoxPlot( Data ):
  Data = { k:v for k,v in Data.items() if (not v is None) and len(v) }

  n = len(Data)
  a = int( np.ceil( np.sqrt( n ) ) )
  b = int( np.ceil( n/a ) )
  fig , axs = plt.subplots( b , a , sharey=True )
  for x in axs.flat[n:]: x.set_visible( False )

  # Fill the plots 
  for (Classification , (Diseases , lData) ) , ax1 in tqdm.tqdm( list( zip( Data.items() , fig.axes ) ) , ncols=Ncol , desc=f"Drawing Box plots" ):
    lStats = GdcStatistics( *lData ) 
    ax1.set_ylim( 1 , 100 )
    box1 = ax1.boxplot( lData , labels= ["$Mutant$" , "$Wild-type$"] , widths= 0.8 , whis=False , showfliers=False , showmeans=True , meanprops=dict(color="grey"), meanline=True, medianprops=dict(color="black") )    
    for i in range( 2 ): ax1.scatter( np.random.normal( i+1 , 0.05 , len( lData[i] ) ) , lData[i] , color=[ "r" , "b" ][i] , alpha=0.5 , s=1 )
    # ax1.text( 0.6 , 30 , f'{Classification}\n$N_{{Mutant}} = {len(lData[0])}$\n$N_{{Wild-type}} = {len(lData[1])}$\n$p_{{value}}={lStats.pvalue:.2e}$' , fontsize="x-small" )

    labels = [] 
    for k,v in sorted( Diseases.items() ):
      k = f"{k} [{v['Mut']}mut|{v['WT']}wt]"
      if a > 3 and len( k ) > 40 :
        index = k.rfind( ' ' , 0 , 40 )
        k = k[:index] + '\n' + k[index:]
      labels.append( k )     

    ax1.text( 0.6 , 25 , "\n".join( sorted( labels ) ) + f'\n$p_{{value}}={lStats.pvalue:.2e}$' , fontsize="x-small" )

  #Add the common y-axis label
  plt.yscale( "log" )
  fig.add_subplot(111, frameon=False)
  plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
  plt.ylabel( "DRG2 TPM-unstranded" , style='italic' )
    
  # Draw the images
  fig.set_size_inches( 16 , 20 )
  plt.tight_layout()
  fig.subplots_adjust( hspace = 0.2 , wspace = 0.0 )
  plt.savefig( args.dest )  
# ======================================================================================================


# ====================================================================================================== 
if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
CacheFile = f".cache/BoxPlot.{args.mutation}.BroadClasses.pkl.gz"
LoadAndClassify( args.src , Classifier , BoxPlot_ForEachClass , DrawBoxPlot , cachefile=CacheFile )    
# ======================================================================================================
