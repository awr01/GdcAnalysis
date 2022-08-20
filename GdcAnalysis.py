from GdcLib import *
from GdcUtils import *
from scipy.stats import ttest_ind
from numpy import mean , var , log10 , warnings
from openpyxl import Workbook
warnings.filterwarnings('ignore')


# ======================================================================================================
def Analyze( aCase ):    
  index = 1                                                                      # Indices chosen for consistency with old code
  if MutationOfInterest in aCase.Mutations: 
    if aCase.Mutations[ MutationOfInterest ].Classification == SilentOrSplice : return
    index = 0

  lDiseaseType = str( aCase.DiseaseType )
  if not lDiseaseType in lCases2: lCases2[ lDiseaseType ] = [ [] , [] ]
  lCases2[ lDiseaseType ][ index ].append( aCase )
# ======================================================================================================

# # ======================================================================================================
# def FinallyTsv():
  # with open( args.dest , "w" ) as dest: 
    # dest.write( f"Disease-type\tGene\tGene-type\tFlagged\tMut-count\tMut-mean\tMut-var\tWT-count\tWT-mean\tWT-var\tMut mean/WT mean\tlog_1.5(ratio)\tt-score\tp-value\tneg log_10(p-value)\n" ) # Write headers
    
    # for DiseaseType , GeneList in Genes.items():  
      # if len( GeneList[0][0] ) == 0 : continue
      # for GeneName , GeneData in tqdm.tqdm( zip( StarCounts.GeneCatalogue , GeneList ) , ncols=Ncol , total=len(GeneList) , desc="Saving" ):
        # num0 , mean0 , var0 = len( GeneData[0] ) , mean( GeneData[0] ) , var( GeneData[0] )               # Calculate the mean and variance of the muts
        # num1 , mean1 , var1 = len( GeneData[1] ) , mean( GeneData[1] ) , var( GeneData[1] )               # Calculate the mean and variance of the WTs
        # if mean0 == 0 or mean1 == 0 : continue

        # meanratio = mean0 / mean1
        # logmeanratio = log10( meanratio ) / log10( 1.5 ) 
        
        # flag = "*" if GeneName in args.genes else ""
        # t01, p01 = ttest_ind( GeneData[0] , GeneData[1] , equal_var = False ) # Calculate the t-score between each pair of lists

        # GeneType = StarCounts.GeneCatalogue[ GeneName ].type
        # dest.write( f"{DiseaseType}\t{GeneName}\t{GeneType}\t{flag}\t{num0}\t{mean0}\t{var0}\t{num1}\t{mean1}\t{var1}\t{meanratio}\t{logmeanratio}\t{t01}\t{p01}\t{-log10(p01)}\n" )
# # ======================================================================================================

# ======================================================================================================
def FinallyXlsx():

  wb = Workbook()
  wb.remove_sheet( wb.active ) # Delete default sheet

  for lDiseaseType , lCases in tqdm.tqdm( sorted( lCases2.items() ) , ncols=Ncol , desc="Analysing" ):
    lMut , lWt = lCases[0] , lCases[1]
    if len( lMut ) == 0 or len( lWt ) == 0 : continue

    ws = wb.create_sheet( lDiseaseType )
    ws.append( [ "Gene" , "Gene-type" , "Mut-count" , "Mut-mean" , "Mut-var" , "WT-count" , "WT-mean" , "WT-var" , "Mut mean/WT mean" , "log_1.5(ratio)" , "t-score" , "p-value" , "neg log_10(p-value)\n" ] ) # Write headers
        
    for GeneName , Gene in tqdm.tqdm( StarCounts.GeneCatalogue.items() , leave=False , ncols=Ncol , desc = lDiseaseType ):
      if Gene.type != "protein_coding" : continue

      lMut2 = []      
      for j in lMut:
        for i in j.StarCounts: 
          lMut2.append( i.TpmUnstranded[ Gene.index ] )
      num0 , mean0 , var0 = len( lMut2 ) , mean( lMut2 ) , var( lMut2 )               # Calculate the mean and variance of the muts
      if mean0 == 0 : continue
      
      lWt2 = []
      for j in lWt:
        for i in j.StarCounts: 
          lWt2.append( i.TpmUnstranded[ Gene.index ] )
      num1 , mean1 , var1 = len( lWt2 )  , mean( lWt2 )  , var( lWt2 )               # Calculate the mean and variance of the WTs
      if mean1 == 0 : continue

      meanratio = mean0 / mean1     
      t01, p01 = ttest_ind( lMut2 , lWt2 , equal_var = False ) # Calculate the t-score between each pair of lists

      ws.append( [ GeneName , Gene.type , num0 , mean0 , var0 , num1 , mean1 , var1 , meanratio , log10( meanratio ) / log10( 1.5 ) , t01 , p01 , -log10(p01) ] )

  wb.save( args.dest )      
# ======================================================================================================



# ======================================================================================================
if len( args.mutations ) != 1 : raise Exception( "Exactly 1 mutation must be specified on the commandline" )
MutationOfInterest , lCases2 = args.mutations[0] , {}         

with open( args.dest , "w" ) as dest: pass



# if args.dest.endswith( ".tsv" ):
# LoadAndForEach( args.src , Analyze , After = FinallyTsv ) # Load src and analyze on-the-fly, making use of the optional After function
# else:
LoadAndForEach( args.src , Analyze , After = FinallyXlsx ) # Load src and analyze on-the-fly, making use of the optional After function
# ======================================================================================================




  