from GdcLib import *
from GdcUtils import *
from openpyxl import Workbook
import matplotlib.pyplot as plt
import numpy as np

# ======================================================================================================
def SortByDisease( aCase ):    
  index = 1                                                                      # Indices chosen for consistency with old code
  if MutationOfInterest in aCase.Mutations: 
    if aCase.Mutations[ MutationOfInterest ].Classification == SilentOrSplice : return
    index = 0

  lDiseaseType = str( aCase.DiseaseType )
  if not lDiseaseType in lCasesByDisease: lCasesByDisease[ lDiseaseType ] = [ [] , [] ]
  lCasesByDisease[ lDiseaseType ][ index ].append( aCase )
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
def SplitString( aStr ):
  lRet , cnt = "" , 0
  for i in aStr:
    if cnt > 23: 
      if i == " " : 
        lRet += "\n"
        cnt = 0
      else:
        lRet += i
        cnt += 1        
    else:
      lRet += i
      cnt += 1
      
  return lRet  
# ======================================================================================================
  

# ======================================================================================================
def Finally():

  DRG2 = StarCounts.GeneCatalogue[ "DRG2" ].index    

  # Filter out empty disease-types
  lCasesByDisease2 = {}
  for lDiseaseType , lCases in lCasesByDisease.items():
    lMut , lWt =  Flatten( lCases[0] , DRG2 ) , Flatten( lCases[1] , DRG2 ) 
    if len( lMut ) > 5 and len( lWt ) > 5 : lCasesByDisease2[ lDiseaseType ] = ( [ lMut , lWt ] , GdcStatistics( lMut , lWt ) )

  # Create the canvas
  fig = plt.gcf()
  ncols = 6
  nrows = int( np.ceil( len( lCasesByDisease2 ) / ncols ) )
  fig , axs = plt.subplots( nrows , ncols , sharey=True )
  plt.yscale( "log" )

  # Fill the plots
  lCasesByDisease2 = { i : j for i,j in sorted( lCasesByDisease2.items(), key=lambda Iter : Iter[1][1].pvalue ) }
  
  for lDiseaseIt , ax1 in zip( lCasesByDisease2.items() , fig.axes ):
    lDiseaseType , lDataStats = lDiseaseIt
    lData , lStats = lDataStats
    
    ax1.set_xlabel( SplitString( lDiseaseType ) )
    ax1.set_ylim( 0 , 100 )
    box1 = ax1.boxplot( lData , labels= ["Mutant" , "Wild-type"] , widths= 0.8 , whis=False , showfliers=False , showmeans=True , meanprops=dict(color="black"), meanline=True, medianprops=dict(visible=False) )    
    for i in range( 2 ): ax1.scatter( np.random.normal( i+1 , 0.05 , len( lData[i] ) ) , lData[i] , color=[ "r" , "b" ][i] , alpha=0.5 , s=1 )

    ax1.text( 0.6 , 50 , f'$N_{{Mutant}} = {len(lData[0])}$\n$N_{{Wild-type}} = {len(lData[1])}$\n$p_{{value}}={lStats.pvalue:.2e}$' , fontsize="x-small" )

  for ax1 in fig.axes[ len( lCasesByDisease2 ): ] : ax1.set_axis_off()


  fig.add_subplot(111, frameon=False)
  plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
  plt.ylabel("TPM-unstranded")
    
  # Draw the images
  fig.set_size_inches( 16 , 4*nrows )
  plt.tight_layout()
  fig.subplots_adjust( hspace = 0.2 , wspace = 0.0 )
  plt.savefig( 'DRG2.pdf' )  
# ======================================================================================================


# ======================================================================================================
if len( args.mutations ) != 1 : raise Exception( "Exactly 1 mutation must be specified on the commandline" )
MutationOfInterest , lCasesByDisease = args.mutations[0] , {}         
LoadAndForEach( args.src , SortByDisease , After = Finally ) # Load src and analyze on-the-fly, making use of the optional After function
# ======================================================================================================
  