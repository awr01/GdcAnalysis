from GdcUtils import *
from GdcPlotUtils import *

# ====================================================================================================== 
if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
GdcAnalysis( lambda aCase: str( aCase.DiseaseType ) , "PerDisease" , DrawVolcanos )
# ======================================================================================================
  