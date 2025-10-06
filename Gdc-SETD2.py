from GdcLib import *
from GdcUtils import *
from openpyxl import Workbook
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import tqdm, argparse

matplotlib.rcParams.update({'font.size': 5})

# ======================================================================================================
def Common_ForEachClass( Class , Cases , index ):
  lMut , lWt , Diseases = SeparateMutantAndWildType( Cases , MutationOfInterest)
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

  for Classification , (_ , lData) in tqdm.tqdm( Data.items() , ncols=Ncol , desc=f"Exporting Xlsx" ):

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

  for (_ , (Diseases , lData) ) , ax1 in tqdm.tqdm( list( zip( Data.items() , fig.axes ) ) , ncols=Ncol , desc=f"Drawing Volcanos" ):

    x0 , y0 , x1 , y1 , x2 , y2 = [] , [] , [] , [] , [] , []

    for GeneName, Stats in lData.items():
      if ( GeneName == "TERT" ):  
        x0.append( Stats.log_mean_ratio_with_error[0] )
        y0.append( Stats.neg_log_pvalue )   
      elif ( Stats.neg_log_pvalue < 2.3 ) or ( np.fabs( Stats.log_mean_ratio_with_error[0] ) < 1 ):  
        x1.append( Stats.log_mean_ratio_with_error[0] )
        y1.append( Stats.neg_log_pvalue )
      else:
        x2.append( Stats.log_mean_ratio_with_error[0] )
        y2.append( Stats.neg_log_pvalue )   

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
  Data = SeparateMutantionType( Cases , MutationOfInterest )

  TERT = StarCounts.GeneCatalogue[ "TERT" ].index

  for k,v in Data.items():
    for K,V in v.items():
      Data[k][K] = FlattenTpmUnstranded( V , TERT )

  return Data
# ======================================================================================================

# ======================================================================================================
def DrawBoxPlot( Data ):

  n = len(Data)
  a = min( 3 , int( np.ceil( np.sqrt( n ) ) ) )
  b = int( np.ceil( n/a ) )
  fig , axs = plt.subplots( b , a , sharey=True )
  if n!=1 :
    for x in axs.flat[n:]: x.set_visible( False )


  # Fill the plots 
  for ( Class , Diseases ) , ax1 in tqdm.tqdm( list( zip( Data.items() , fig.axes ) ) , ncols=Ncol , desc=f"Drawing Box plots" ):

    D = {}
    for k,v in Diseases.items():
      for K,V in v.items():
        if not K in D: D[K] = []
        D[K].extend( V )

    ax1.set_ylim( 1e-2 , 1e3 )

    keys = [ f"{k}\n[{len(v)}]" for k,v in sorted( D.items() ) ]
    vals = [ v for k,v in sorted( D.items() ) ]

    box1 = ax1.boxplot( vals , tick_labels=keys , widths= 0.8 , whis=False , showfliers=False , showmeans=True , meanprops=dict(color="grey"), meanline=True, medianprops=dict(color="black") )    

    for i,v in enumerate(vals): 
      ax1.scatter( np.random.normal( i+1 , 0.05 , len( v ) ) , v , color=[ "r" , "b" ][i==len(vals)-1] , alpha=0.5 , s=1 )

      if i!=len(vals)-1 : 
        try:
          lStats = GdcStatistics( v , vals[-1] ) 
          ax1.text( i+0.6 , 100 + (100*(i%2)) , f'$p_{{value}}={lStats.pvalue:.2e}$' )
        except: pass

    ax1.text( 0.6 , 300 , "\n".join( sorted( Diseases.keys() ) ) )

  #Add the common y-axis label
  plt.yscale( "log" )
  fig.add_subplot(111, frameon=False)
  plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
  plt.ylabel( "TERT TPM-unstranded" , style='italic' )
    
  # Draw the images
  fig.set_size_inches( 16 , 20 )
  plt.tight_layout()
  fig.subplots_adjust( hspace = 0.2 , wspace = 0.0 )
  plt.savefig( args.dest )  
# ======================================================================================================




# ======================================================================================================
def ScatterPlot_ForEachClass( Class , Cases , index ):
  return FlattenTpmUnstranded( Cases , StarCounts.GeneCatalogue[ "SETD2" ].index ),  FlattenTpmUnstranded( Cases , StarCounts.GeneCatalogue[ "TERT" ].index )
# ======================================================================================================

# ======================================================================================================
def DrawScatterPlot( Data ):

  n = len(Data)
  a = min( 3 , int( np.ceil( np.sqrt( n ) ) ) )
  b = int( np.ceil( n/a ) )
  fig , axs = plt.subplots( b , a , sharey=True )
  if n!=1 : 
    for x in axs.flat[n:]: x.set_visible( False )

  # Fill the plots 
  for ( Class , Data ) , ax1 in tqdm.tqdm( list( zip( Data.items() , fig.axes ) ) , ncols=Ncol , desc=f"Drawing Scatter plots" ):
  #   ax1.set_xlim( 1 , 1e3 )
  #   ax1.set_ylim( 1e-4 , 1e2 )
  #   ax1.scatter( x=Data[0] , y=Data[1] , s=1 )    
  #   ax1.set_xscale( "log" )
  #   ax1.text( 2 , 30 , Class )

  #   z = np.polyfit( Data[0] , Data[1], 1)
  #   p = np.poly1d(z)
  #   ax1.plot( [1,1e3] , p([1,1e3]) , "r--" )

    Data = sorted( zip( *Data ) )
    x = int( np.around( len(Data) / 10 ) )
    Data = [ [y[1] for y in Data[:x] ] , [y[1] for y in Data[x:] ] ]

    ax1.set_ylim( 1e-4 , 1e2 )
    box1 = ax1.boxplot( Data , tick_labels=[ "Bottom 10%" , "Top 90%" ] , widths= 0.8 , whis=False , showfliers=False , showmeans=True , meanprops=dict(color="grey"), meanline=True, medianprops=dict(color="black") )    

    for i,v in enumerate(Data): 
      ax1.scatter( np.random.normal( i+1 , 0.05 , len( v ) ) , v , color=[ "r" , "b" ][i] , alpha=0.5 , s=1 )

    ax1.text( .6 , 45 , Class )
    try:
      lStats = GdcStatistics( *Data ) 
      ax1.text( 0.6 , 20 , f'$p_{{value}}={lStats.pvalue:.2e}$' )
    except: pass

  #Add the common y-axis label
  plt.yscale( "log" )
  fig.add_subplot(111, frameon=False)
  plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
  plt.ylabel( "TERT TPM-unstranded" , style='italic' )
  # plt.xlabel( "SETD2 TPM-unstranded" , style='italic' )

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
MutationOfInterest = "SETD2"

parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball' )
parser.add_argument( '--dest' , help='The destination file' )
parser.add_argument( '--output' , required=True , choices=[ 'Excel' , 'Volcano' , 'BoxPlot' , 'ScatterPlot' ] , help='The output type' )
parser.add_argument( '--classification' , required=True , choices=[ 'PerDisease' , 'BroadClasses' , 'OneBin' ] , help='Treat per-pisease or use broader classes' )

args = parser.parse_args()

if not args.src .endswith( ".tar"  ): raise Exception( "Source file must have '.tar' file-extension" )


if args.dest is None : 
  if    args.output == "Excel": args.dest = f"SETD2-TERT-{args.classification}.xlsx"
  else:                         args.dest = f"SETD2-TERT-{args.output}-{args.classification}.pdf"
  print( f"Set destination to '{args.dest}" )

if args.classification == 'PerDisease':
  classifierfn , maxthreads = lambda aCase: str( aCase.DiseaseType ) , None
elif args.classification == 'BroadClasses':                                 
  classifierfn , maxthreads = BroadClasses , 3
elif args.classification == 'OneBin':                                 
  classifierfn , maxthreads = lambda aCase: "Everything" , 1

if args.output == "Excel":
  if not args.dest.endswith( ".xlsx" ): raise Exception( "Destination file must have '.xlsx' file-extension" )
  cacheprefix , foreachfn , exportfn = "Common" , Common_ForEachClass , ExportXlsx
elif args.output == "Volcano":
  if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
  cacheprefix , foreachfn , exportfn = "Common" , Common_ForEachClass , DrawVolcanos
elif args.output == "BoxPlot": # BoxPlot
  if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
  cacheprefix , foreachfn , exportfn = "BoxPlot" , BoxPlot_ForEachClass , DrawBoxPlot
else: #ScatterPlot
  if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
  cacheprefix , foreachfn , exportfn = "Scatter" , ScatterPlot_ForEachClass , DrawScatterPlot


if not os.path.isdir( ".cache" ): os.mkdir( ".cache" )
CacheFile = f".cache/SETD2-TERT-{cacheprefix}-{args.classification}.pkl.gz"
LoadAndClassify( args.src , classifierfn , foreachfn , exportfn , cachefile=CacheFile , maxthreads=maxthreads )    
# ======================================================================================================
