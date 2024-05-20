from GdcLib import *
from GdcUtils import *
from openpyxl import Workbook
import matplotlib.pyplot as plt
import numpy as np
import tqdm, argparse

# ======================================================================================================
def Common_ForEachClass( Class , Cases , index ):
  lMut , lWt , Diseases = SeparateMutandAndWildType( Cases , MutationOfInterest)
  if len( lMut ) == 0 or len( lWt ) == 0 : return
  
  Results = {}
  for GeneName , Gene in tqdm.tqdm( sorted( StarCounts.GeneCatalogue.items() ) , leave=False , ncols=Ncol , desc=f"{Class}: Analysing", position=index ):
    if Gene.type != "protein_coding" : continue
    lRet = GdcStatistics( FlattenTpmUnstranded( lMut , Gene.index ) , FlattenTpmUnstranded( lWt , Gene.index ) )
    if lRet is None: continue
    if np.isnan( lRet.neg_log_pvalue ) : continue
    Results[ GeneName ] = lRet
  return Diseases , Results
# ======================================================================================================

# ======================================================================================================
def ExportXlsx( Data ): 
  wb , ws = Workbook() , None
  wb.remove_sheet( wb.active ) # Delete default sheet

  Data = { k:v for k,v in Data.items() if (not v is None) and len(v) }

  for Classification , (Diseases , lData) in tqdm.tqdm( Data.items() , ncols=Ncol , desc=f"Exporting Xlsx" ):
    
    ws = wb.create_sheet( Classification )
    ws.append( [ "Gene" , 
                  "Mut-count" , "Mut-mean" , "Mut-mean-error" , "Mut-std.dev" , 
                  "WT-count"  , "WT-mean"  , "WT-mean-error"  , "WT-std.dev" , 
                  "Mut mean/WT mean" , "error(Mut mean/WT mean)" , 
                  "log_1.5(ratio)" , "error(log_1.5(ratio))" , 
                  "t-score" , "p-value" , "neg log_10(p-value)" ] ) # Write headers

    for GeneName, aStats in sorted( lData.items() , key = lambda x : x[1].pvalue ):
      ws.append( [ GeneName , 
                    aStats.mut.count , aStats.mut.mean , aStats.mut.mean_error , aStats.mut.sd , 
                    aStats.wt.count  , aStats.wt.mean  , aStats.wt.mean_error  , aStats.wt.sd , 
                    aStats.mean_ratio_with_error[0] , aStats.mean_ratio_with_error[1] , 
                    aStats.log_mean_ratio_with_error[0] , aStats.log_mean_ratio_with_error[1] , 
                    aStats.tscore , aStats.pvalue , aStats.neg_log_pvalue ] )

  wb.save( args.dest ) 
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
def BoxPlot_ForEachClass( Class , Cases , index ):
  lMut , lWt , Diseases = SeparateMutandAndWildType( Cases , MutationOfInterest )
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
def BroadClasses( aCase ):    
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

MutationOfInterest = "ATRX"

parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball' )
parser.add_argument( '--dest' , help='The destination file' )
# parser.add_argument( '--mutation' , required=True , help='Mutation of interest')
parser.add_argument( '--output' , required=True , choices=[ 'Excel' , 'Volcano' , 'BoxPlot' ] , help='The output type' )
parser.add_argument( '--classification' , required=True , choices=[ 'PerDisease' , 'BroadClasses' ] , help='Treat per-pisease or use broader classes' )

args = parser.parse_args()

if not args.src .endswith( ".tar"  ): raise Exception( "Source file must have '.tar' file-extension" )


if args.dest is None : 
  if    args.output == "Excel": args.dest = f"ATRX-DRG2-{args.classification}.xlsx"
  else:                         args.dest = f"ATRX-DRG2-{args.output}-{args.classification}.pdf"
  print( f"Set destination to '{args.dest}" )

if args.classification == 'PerDisease':
  classifierfn , maxthreads = lambda aCase: str( aCase.DiseaseType ) , None
else: # BroadClasses                                 
  classifierfn , maxthreads = BroadClasses , 3

if args.output == "Excel":
  if not args.dest.endswith( ".xlsx" ): raise Exception( "Destination file must have '.xlsx' file-extension" )
  cacheprefix , foreachfn , exportfn = "Common" , Common_ForEachClass , ExportXlsx
elif args.output == "Volcano":
  if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
  cacheprefix , foreachfn , exportfn = "Common" , Common_ForEachClass , DrawVolcanos
else: # BoxPlot
  if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
  cacheprefix , foreachfn , exportfn = "BoxPlot" , BoxPlot_ForEachClass , DrawBoxPlot

if not os.path.isdir( ".cache" ): os.mkdir( ".cache" )
CacheFile = f".cache/ATRX-DRG2-{cacheprefix}-{args.classification}.pkl.gz"
LoadAndClassify( args.src , classifierfn , foreachfn , exportfn , cachefile=CacheFile , maxthreads=maxthreads )    

# ======================================================================================================
