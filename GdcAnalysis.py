# from GdcUtils import *
# from openpyxl import Workbook
# import tqdm

# # ======================================================================================================
# def ExportXlsx( Data ): 
#   wb , ws = Workbook() , None
#   wb.remove_sheet( wb.active ) # Delete default sheet

#   for lDiseaseType , lData in tqdm.tqdm( Data.items() , ncols=Ncol , desc=f"Exporting Xlsx" ):

#     if lData is None: continue

#     ws = wb.create_sheet( lDiseaseType )
#     ws.append( [ "Gene" , 
#                   "Mut-count" , "Mut-mean" , "Mut-mean-error" , "Mut-std.dev" , 
#                   "WT-count"  , "WT-mean"  , "WT-mean-error"  , "WT-std.dev" , 
#                   "Mut mean/WT mean" , "error(Mut mean/WT mean)" , 
#                   "log_1.5(ratio)" , "error(log_1.5(ratio))" , 
#                   "t-score" , "p-value" , "neg log_10(p-value)" ] ) # Write headers

#     for GeneName, aStats in sorted( lData.items() , key = lambda x : x[1].pvalue ):
#       ws.append( [ GeneName , 
#                     aStats.mut.count , aStats.mut.mean , aStats.mut.mean_error , aStats.mut.sd , 
#                     aStats.wt.count  , aStats.wt.mean  , aStats.wt.mean_error  , aStats.wt.sd , 
#                     aStats.mean_ratio_with_error[0] , aStats.mean_ratio_with_error[1] , 
#                     aStats.log_mean_ratio_with_error[0] , aStats.log_mean_ratio_with_error[1] , 
#                     aStats.tscore , aStats.pvalue , aStats.neg_log_pvalue ] )

#   wb.save( args.dest ) 
# # ======================================================================================================

# # ====================================================================================================== 
# if not args.dest.endswith( ".xlsx" ): raise Exception( "Destination file must have '.xlsx' file-extension" )
# GdcAnalysis( lambda aCase: str( aCase.DiseaseType ) , "PerDisease" , ExportXlsx )
# # ======================================================================================================
