# from GdcLib import *
from GdcUtils import *
from openpyxl import Workbook
# import numpy as np
# import os , hashlib , bz2 , _pickle

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
  
#   Results = {}
#   for GeneName , Gene in tqdm.tqdm( sorted( StarCounts.GeneCatalogue.items() ) , leave=False , ncols=Ncol , desc=f"{Class}: Analysing", position=index ):
#     if Gene.type != "protein_coding" : continue
#     lRet = GdcStatistics( Flatten( lMut , Gene.index ) , Flatten( lWt , Gene.index ) )
#     if lRet is None: continue
#     if np.isnan( lRet.neg_log_pvalue ) : continue

#     Results[ GeneName ] = lRet

#   return Results
# # ======================================================================================================


# ======================================================================================================
def ExportXlsx( Data ): 
  wb , ws = Workbook() , None
  wb.remove_sheet( wb.active ) # Delete default sheet

  for lDiseaseType , lData in tqdm.tqdm( Data.items() , ncols=Ncol , desc=f"Exporting Xlsx" ):

    if lData is None: continue

    ws = wb.create_sheet( lDiseaseType )
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
# LoadAndClassify( args.src , lambda aCase: str( aCase.DiseaseType ) , ForEachClass , ExportXlsx )
# ======================================================================================================

GdcAnalysis( lambda aCase: str( aCase.DiseaseType ) , "PerDisease" , ExportXlsx )