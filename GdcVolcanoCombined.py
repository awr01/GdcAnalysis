from GdcUtils import *
from GdcPlotUtils import *

# ======================================================================================================
def Classifier( aCase ):    
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
if not args.dest.endswith( ".pdf" ): raise Exception( "Destination file must have '.pdf' file-extension" )
GdcAnalysis( Classifier , "BroadClasses" , DrawVolcanos , maxthreads=3 )
# ======================================================================================================  