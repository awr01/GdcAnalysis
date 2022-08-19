from GdcLib import *
from GdcUtils import *
from scipy.stats import ttest_ind
from numpy import mean , var , log10 , warnings
warnings.filterwarnings('ignore')


from openpyxl import Workbook
wb = Workbook()

wb.remove_sheet( wb.active ) # Delete default sheet

# ======================================================================================================
def Analyze( aCase ):    

  if not aCase.DiseaseType in Genes:
    Genes[ aCase.DiseaseType ] = [ ( [],[] ) for i in StarCounts.GeneCatalogue ]

  index = 1                                                                      # Indices chosen for consistency with old code
  if MutationOfInterest in aCase.Mutations: 
    if aCase.Mutations[ MutationOfInterest ].Classification == SilentOrSplice : return       # Ignore silents and slices
    index = 0
        
  for lStarCount in aCase.StarCounts:                                            # Then iterate over each star-count file
    for TpmUnstranded , GeneData in zip( lStarCount.TpmUnstranded , Genes[ aCase.DiseaseType ] ):
      if not TpmUnstranded is None: GeneData[ index ].append( TpmUnstranded ) 
# ======================================================================================================

# ======================================================================================================
def FinallyTsv():
  with open( args.dest , "w" ) as dest: 
    for DiseaseType , GeneList in Genes.items():
      dest.write( f"Disease-type\tGene\tGene-type\tFlagged\tMut-count\tMut-mean\tMut-var\tWT-count\tWT-mean\tWT-var\tMut mean/WT mean\tlog_1.5(ratio)\tt-score\tp-value\tneg log_10(p-value)\n" ) # Write headers
      
      for GeneName , GeneData in tqdm.tqdm( zip( StarCounts.GeneCatalogue , GeneList ) , ncols=Ncol , total=len(GeneList) , desc="Saving" ):
        num0 , mean0 , var0 = len( GeneData[0] ) , mean( GeneData[0] ) , var( GeneData[0] )               # Calculate the mean and variance of the muts
        num1 , mean1 , var1 = len( GeneData[1] ) , mean( GeneData[1] ) , var( GeneData[1] )               # Calculate the mean and variance of the WTs
        if mean0 == 0 or mean1 == 0 : continue

        meanratio = mean0 / mean1
        logmeanratio = log10( meanratio ) / log10( 1.5 ) 
        
        flag = "*" if GeneName in args.genes else ""
        t01, p01 = ttest_ind( GeneData[0] , GeneData[1] , equal_var = False ) # Calculate the t-score between each pair of lists

        GeneType = StarCounts.GeneCatalogue[ GeneName ].type
        dest.write( f"{DiseaseType}\t{GeneName}\t{GeneType}\t{flag}\t{num0}\t{mean0}\t{var0}\t{num1}\t{mean1}\t{var1}\t{meanratio}\t{logmeanratio}\t{t01}\t{p01}\t{-log10(p01)}\n" )
# ======================================================================================================

# ======================================================================================================
def FinallyXlsx():
  for DiseaseType , GeneList in Genes.items():
    ws = wb.create_sheet( DiseaseType )
    ws.append( [ "Gene" , "Gene-type" , "Flagged" , "Mut-count" , "Mut-mean" , "Mut-var" , "WT-count" , "WT-mean" , "WT-var" , "Mut mean/WT mean" , "log_1.5(ratio)" , "t-score" , "p-value" , "neg log_10(p-value)\n" ] ) # Write headers
    
    for GeneName , GeneData in tqdm.tqdm( zip( StarCounts.GeneCatalogue , GeneList ) , ncols=Ncol , total=len(GeneList) , desc="Saving" ):
      num0 , mean0 , var0 = len( GeneData[0] ) , mean( GeneData[0] ) , var( GeneData[0] )               # Calculate the mean and variance of the muts
      num1 , mean1 , var1 = len( GeneData[1] ) , mean( GeneData[1] ) , var( GeneData[1] )               # Calculate the mean and variance of the WTs
      if mean0 == 0 or mean1 == 0 : continue

      meanratio = mean0 / mean1
      logmeanratio = log10( meanratio ) / log10( 1.5 ) 
      
      flag = "*" if GeneName in args.genes else ""
      t01, p01 = ttest_ind( GeneData[0] , GeneData[1] , equal_var = False ) # Calculate the t-score between each pair of lists

      GeneType = StarCounts.GeneCatalogue[ GeneName ].type
      ws.append( [ GeneName , GeneType , flag , num0 , mean0 , var0 , num1 , mean1 , var1 , meanratio , logmeanratio , t01 , p01 , -log10(p01) ] )
  wb.save( args.dest )      
# ======================================================================================================



# ======================================================================================================
if len( args.mutations ) != 1 : raise Exception( "Exactly 1 mutation must be specified on the commandline" )
MutationOfInterest , Genes = args.mutations[0] , {}         

with open( args.dest , "w" ) as dest: pass

if args.dest.endswith( ".tsv" ):
  LoadAndForEach( args.src , Analyze , After = FinallyTsv ) # Load src and analyze on-the-fly, making use of the optional After function
else:
  LoadAndForEach( args.src , Analyze , After = FinallyXlsx ) # Load src and analyze on-the-fly, making use of the optional After function



# ======================================================================================================
