from GdcLib import *
from GdcPlotUtils import *

# ====================================================================================================== 
if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
CacheFile = f".cache/BoxPlot.{args.mutation}.PerDisease.pkl.gz"
LoadAndClassify( args.src , lambda aCase: str( aCase.DiseaseType ) , BoxPlot_ForEachClass , DrawBoxPlot , cachefile=CacheFile )    
# ======================================================================================================
