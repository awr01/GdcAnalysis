import argparse , requests, json, tarfile, io, codecs, gzip , tqdm , os , hashlib , bz2 , _pickle
from GdcLib import *

utf8reader = codecs.getreader( 'utf-8' )


# Read the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--dest' , required=True , help='The destination tgz file')
args = parser.parse_args()
if not args.dest.endswith( ".tgz" ): raise Exception( "Destination file must have '.tgz' file-extension" )


if not os.path.isdir( ".cache" ): os.mkdir( ".cache" )

# ======================================================================================================
# Magic to simplify writing filters
def And( *args ): return { "op":"and" , "content":[ i for i in args ] }
def Or( *args ):  return { "op":"or" , "content":[ i for i in args ] }
def Eq( key , value ):  return { "op":"=" , "content":{ "field" : key , "value": value } }
def In( key , value ):  return { "op":"in" , "content":{ "field" : key , "value": value } }
# ======================================================================================================


# ======================================================================================================
def CachedGetJson( aParams ):
  lParams = hashlib.md5( str(aParams).encode('utf-8') ).hexdigest()
  lFile = f".cache/{lParams}.pkl.bz"
  
  if os.path.isfile( lFile ):
    print( "> Using cached query" , flush=True )
    with bz2.BZ2File( lFile , 'r' ) as src: return _pickle.load( src )  

  print( "> Downloading query" , flush=True )
  lResponse = requests.get( "https://api.gdc.cancer.gov/files" , params = aParams )  
  if lResponse.status_code != 200: raise Exception( f"Requests returned status-code {lResponse.status_code}" )  
  lResponse = json.loads( lResponse.content.decode("utf-8") ) # Get response and convert to python
  lResponse = lResponse["data"]["hits"]
    
  with bz2.BZ2File( lFile , 'w' ) as dest: _pickle.dump( lResponse , dest )

  return lResponse
# ======================================================================================================

# ======================================================================================================
def GetJson( aFilter , aFields = [] ):
  return CachedGetJson( { "filters":json.dumps( aFilter ) , "fields":( ",".join(aFields) ) , "size": "99999999" } )
# ======================================================================================================



# ======================================================================================================
lInfo = [ "cases.case_id" , "cases.diagnoses.age_at_diagnosis" , "cases.project.disease_type" , "experimental_strategy" , "access" ]       
lFileInfo = GetJson( Or( And( Eq( "files.experimental_strategy" , "RNA-Seq" ), Eq( "files.analysis.workflow_type" , "STAR - Counts"), Eq( "files.data_type" , "Gene Expression Quantification") ) ,       
                         And( Eq( "files.experimental_strategy" , "WXS" ) ) ) , lInfo )

print( "#Files = " , len( lFileInfo ) )

OpenBoth = []
OpenRnaSeq = []
OpenWXS = []
Controlled = []



lCases = {}

for j in tqdm.tqdm( lFileInfo , ncols=Ncol , desc="Collating available Data" ):  
  CaseId = j["cases"][0]["case_id"]
  FileId = j["id"]
  ExperimentalStrategy = j["experimental_strategy"]
  Access = j["access"]

  if not CaseId in lCases: lCases[ CaseId ] = {}
  if not Access in lCases[ CaseId ]: lCases[ CaseId ][ Access ] = {}
  if not ExperimentalStrategy in lCases[ CaseId ][ Access ] : lCases[ CaseId ][ Access ][ ExperimentalStrategy ] = []

  lCases[ CaseId ][ Access ][ ExperimentalStrategy ].append( FileId )



for CaseId , Access in lCases.items():

  if "open" in Access:
    Types = Access[ "open" ]
    if len( Types ) == 2 : OpenBoth.append( CaseId )
    elif "WXS" in Types : OpenWXS.append( CaseId )
    else: OpenRnaSeq.append( CaseId )
  else:
    Controlled.append( CaseId )


print( f"Cases with both WXS and RnaSeq data open:                        { len(OpenBoth) }" )
print( f"Cases with WXS data open (RnaSeq controlled or absent):          { len(OpenWXS) }" )
print( f"Cases with RnaSeq data open (WXS controlled or absent):          { len(OpenRnaSeq)}" )
print( f"Cases with only controlled data (WXS only, RnaSeq only or both): { len(Controlled)}" )
