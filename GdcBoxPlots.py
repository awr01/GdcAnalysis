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


DiseaseTypeLut = {
  "Brain Lower Grade Glioma" : "Brain Lower\nGrade Glioma" ,
  "Sarcoma" : "Sarcoma" ,
  "Glioblastoma Multiforme" : "Glioblastoma Multiforme" ,
  "Pheochromocytoma and Paraganglioma" : "Pheochromocytoma\n& Paraganglioma" ,
  "Neuroblastoma" : "Neuroblastoma" ,
  "Uterine Corpus Endometrial Carcinoma" : "Uterine Corpus\nEndometrial Carcinoma" ,
  "Ovarian Serous Cystadenocarcinoma" : "Ovarian Serous\nCystadenocarcinoma" ,
  "Breast Invasive Carcinoma" : "Breast Invasive Carcinoma" ,
  "Rectum Adenocarcinoma" : "Rectum Adenocarcinoma" ,
  "Kidney Renal Clear Cell Carcinoma" : "Kidney Renal\nClear Cell Carcinoma" ,
  "Bladder Urothelial Carcinoma" : "Bladder Urothelial Carcinoma" ,
  "Skin Cutaneous Melanoma" : "Skin Cutaneous Melanoma" ,
  "Colon Adenocarcinoma" : "Colon Adenocarcinoma" ,
  "Stomach Adenocarcinoma" : "Stomach Adenocarcinoma" ,
  "Lung Squamous Cell Carcinoma" : "Lung Squamous Cell Carcinoma" ,
  "Head and Neck Squamous Cell Carcinoma" : "Head and Neck\nSquamous Cell Carcinoma" ,
  "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma" : "Cervical Squamous Cell Carcinoma\n& Endocervical Adenocarcinoma" ,
  "Lung Adenocarcinoma" : "Lung Adenocarcinoma"
}


# ======================================================================================================
def Finally():

  DRG2 = StarCounts.GeneCatalogue[ "DRG2" ].index    

  # Filter out empty disease-types
  lCasesByDisease2 = {}
  for lDiseaseType , lCases in lCasesByDisease.items():
    lMut , lWt =  Flatten( lCases[0] , DRG2 ) , Flatten( lCases[1] , DRG2 ) 
    lCasesByDisease2[ lDiseaseType ] = ( [ lMut , lWt ] , GdcStatistics( lMut , lWt ) )

  # Create the canvas
  fig = plt.gcf()
  ncols = 6
  nrows = int( np.ceil( len( DiseaseTypeLut ) / ncols ) )
  fig , axs = plt.subplots( nrows , ncols , sharey=True )
  plt.yscale( "log" )

  # Fill the plots 
  for lDiseaseType , ax1 in zip( DiseaseTypeLut , fig.axes ):
    lData , lStats = lCasesByDisease2[ lDiseaseType ] 
    ax1.set_xlabel( DiseaseTypeLut[ lDiseaseType ] , style='italic' , labelpad=1 )
    ax1.set_ylim( 0 , 100 )
    box1 = ax1.boxplot( lData , labels= ["$Mutant$" , "$Wild-type$"] , widths= 0.8 , whis=False , showfliers=False , showmeans=True , meanprops=dict(color="grey"), meanline=True, medianprops=dict(color="black") )    
    for i in range( 2 ): ax1.scatter( np.random.normal( i+1 , 0.05 , len( lData[i] ) ) , lData[i] , color=[ "r" , "b" ][i] , alpha=0.5 , s=1 )
    ax1.text( 0.6 , 50 , f'$N_{{Mutant}} = {len(lData[0])}$\n$N_{{Wild-type}} = {len(lData[1])}$\n$p_{{value}}={lStats.pvalue:.2e}$' , fontsize="x-small" )

  # Blank the remaining subplots
  for ax1 in fig.axes[ len( lCasesByDisease2 ): ] : ax1.set_axis_off()

  #Add the common y-axis label
  fig.add_subplot(111, frameon=False)
  plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
  plt.ylabel( "TPM-unstranded" , style='italic' )
    
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
  