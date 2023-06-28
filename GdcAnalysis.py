from GdcLib import *
from GdcUtils import *
from openpyxl import Workbook


# ======================================================================================================
def SortByDisease( aCase ):    
  index = 1                                                                      # Indices chosen for consistency with old code
  if MutationOfInterest in aCase.Mutations: 
    if aCase.Mutations[ MutationOfInterest ].Classification == SilentOrSplice : return
    index = 0

  lDiseaseType = str( aCase.DiseaseType )
  if not lDiseaseType in lCases2: lCases2[ lDiseaseType ] = [ [] , [] ]
  lCases2[ lDiseaseType ][ index ].append( aCase )
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
  for lDiseaseType , lCases in tqdm.tqdm( sorted( lCases2.items() ) , ncols=Ncol , desc="Analysing" ):
    lMut , lWt = lCases[0] , lCases[1]
    if len( lMut ) == 0 or len( lWt ) == 0 : continue

    HeaderFn( lDiseaseType )
    
    for GeneName , Gene in tqdm.tqdm( sorted( StarCounts.GeneCatalogue.items() ) , leave=False , ncols=Ncol , desc = lDiseaseType ):
      if Gene.type != "protein_coding" : continue

      lRet = GdcStatistics( Flatten( lMut , Gene.index ) , Flatten( lWt , Gene.index ) )
      if lRet is None : continue
      # print( lRet[ "test" , "test.test" ] )

      OutputFn( lDiseaseType , GeneName , Gene.type , GdcStatistics( Flatten( lMut , Gene.index ) , Flatten( lWt , Gene.index ) ))   
# ======================================================================================================

# ======================================================================================================
def HeaderTsv( DiseaseType ): pass

def OutputTsv( DiseaseType , GeneName , GeneType , aStats ):
  dest.write( f"{DiseaseType}\t{GeneName}\t{GeneType}\t" )
  # dest.write( f"{aStats.mut.count}\t{aStats.mut.mean}\t{aStats.mut.mean_error}\t{aStats.mut.sd}\t" )
  # dest.write( f"{aStats.wt.count }\t{aStats.wt.mean }\t{aStats.wt.mean_error }\t{aStats.wt.sd}\t" )
  dest.write( f"{aStats.mean_ratio_with_error[0]}\t{aStats.mean_ratio_with_error[1]}\t")
  dest.write( f"{aStats.log_mean_ratio_with_error[0]}\t{aStats.log_mean_ratio_with_error[1]}\t" )
  dest.write( f"{aStats.tscore}\t{aStats.pvalue}\t{aStats.neg_log_pvalue}\n" )
# ======================================================================================================

# ======================================================================================================
def HeaderXlsx( DiseaseType ):
  global ws
  ws = wb.create_sheet( DiseaseType )
  ws.append( [ "Gene" , "Gene-type" , "Mut-count" , "Mut-mean" , "Mut-mean-error" , "Mut-std.dev" , "WT-count" , "WT-mean" , "WT-mean-error" , "WT-std.dev" , "Mut mean/WT mean" , "log_1.5(ratio)" , "t-score" , "p-value" , "neg log_10(p-value)\n" ] ) # Write headers

def OutputXlsx( DiseaseType , GeneName , GeneType , aStats ):
  if aStats is None: return
  ws.append( [ GeneName , GeneType , 
                # aStats.mut.count , aStats.mut.mean , aStats.mut.mean_error , aStats.mut.sd , 
                # aStats.wt.count  , aStats.wt.mean  , aStats.wt.mean_error  , aStats.wt.sd , 
                aStats.mean_ratio_with_error[0] , aStats.mean_ratio_with_error[1] , 
                aStats.log_mean_ratio_with_error[0] , aStats.log_mean_ratio_with_error[1] , 
                aStats.tscore , aStats.pvalue , aStats.neg_log_pvalue ] )
# ======================================================================================================

# ======================================================================================================
if len( args.mutations ) != 1 : raise Exception( "Exactly 1 mutation must be specified on the commandline" )
MutationOfInterest , lCases2 = args.mutations[0] , {}         

if args.dest.endswith( ".tsv" ):
  with open( args.dest , "w" ) as dest:
    dest.write( f"Disease-type\tGene\tGene-type\t" )
    # dest.write( f"Mut-count\tMut-mean\tMut-mean-error\tMut-std.dev\t" )
    # dest.write( f"WT-count\tWT-mean\tWT-mean-error\tWT-std.dev\t" )
    dest.write( f"Mut mean/WT mean\t" )
    dest.write( f"log_1.5(ratio)\t" )
    dest.write( f"t-score\tp-value\tneg log_10(p-value)\n" ) # Write headers  
    HeaderFn , OutputFn = HeaderTsv , OutputTsv
    LoadAndForEach( args.src , SortByDisease , After = Finally ) # Load src and analyze on-the-fly, making use of the optional After function
else:
  with open( args.dest , "w" ) as dest: pass
  HeaderFn , OutputFn = HeaderXlsx , OutputXlsx
  wb , ws = Workbook() , None
  wb.remove_sheet( wb.active ) # Delete default sheet
  LoadAndForEach( args.src , SortByDisease , After = Finally ) # Load src and analyze on-the-fly, making use of the optional After function
  wb.save( args.dest )    
  
# ======================================================================================================
  