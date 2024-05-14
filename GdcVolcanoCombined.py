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

  return Class
# ======================================================================================================

# ======================================================================================================
def DrawVolcanos( Data ):    
  Data = { k:v for k,v in Data.items() if (not v is None) and len(v) }

  n = len(Data)
  a = int( np.ceil( np.sqrt( n ) ) )
  b = int( np.ceil( n/a ) )
  fig , axs = plt.subplots( b , a , sharey=True )
  for x in axs.flat[n:]: x.set_visible( False )

  for (Classification , (Diseases , lData) ) , ax1 in tqdm.tqdm( list( zip( Data.items() , fig.axes ) ) , ncols=Ncol , desc=f"Drawing Volcanos" ):

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

    # ax1.set_xlabel( Classification , style='italic' , labelpad=1 )
    ax1.set_xlim( -25 , 25 )
    ax1.set_ylim( 1/200 , 200 )
    ax1.scatter( x1 , y1 , color="0.75" , s=1 )
    ax1.scatter( x2 , y2 , color="b" , s=1 )
    ax1.scatter( x0 , y0 , color="r" , s=1 )
    ax1.grid( True )

    labels = [] 
    for k,v in sorted( Diseases.items() ):
      k = f"{k} [{v['Mut']}mut|{v['WT']}wt]"
      if a > 3 and len( k ) > 40 :
        index = k.rfind( ' ' , 0 , 40 )
        k = k[:index] + '\n' + k[index:]
      labels.append( k )     

    ax1.text( -24 , 0.01 , "\n".join( sorted( labels ) ) , fontsize="x-small" )

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
GdcAnalysis( Classifier , "BroadClasses" , DrawVolcanos , maxthreads=3 )
# ======================================================================================================  