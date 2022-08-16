import argparse , requests, json, tarfile, io, codecs, gzip , tqdm , os , hashlib , bz2 , _pickle
from GdcLib import *

utf8reader = codecs.getreader( 'utf-8' )


# Read the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--dest' , required=True , help='The destination tgz file')
args = parser.parse_args()
if not args.dest.endswith( ".tgz" ): raise Exception( "Destination file must have '.tgz' file-extension" )


if not os.path.isdir( ".cache" ): os.path.mkdir( ".cache" )

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
def CachedGetTar( aKeys ):
  lKeys = hashlib.md5( str( sorted( aKeys ) ).encode('utf-8') ).hexdigest()
  lFile = f".cache/{lKeys}.tgz"

  if os.path.isfile( lFile ): 
    print( " Using cached!" , end = "" , flush=True )  
  else:
    print( " Downloading!" , end = "" , flush=True )
    lResponse = requests.post( "https://api.gdc.cancer.gov/data" , data = json.dumps( { "ids":aKeys } ) , headers = { "Content-Type" : "application/json" } )
    if lResponse.status_code != 200: raise Exception( f"Requests returned status-code {lResponse.status_code}" )  
    with open( lFile , "wb" ) as dest: dest.write( lResponse.content )
          
  src = tarfile.open( lFile , mode='r:gz' )
  return src
# ======================================================================================================

# ======================================================================================================
def AddWxsFileToCase( aCase , aFile ):
  for line in utf8reader( gzip.GzipFile( fileobj = aFile ) ):                  # Iterate over each line in the file
    if line[0] == "#" or line.startswith( "Hugo_Symbol" ) : continue           # Ignore comments and headers
    line = [ i.strip() for i in line.split( "\t" , maxsplit = 37 ) ]           # Split the line at tabs up to where we need it      
    if not line[0] in aCase.Mutations: aCase.Mutations[ line[0] ] = Mutation()
    aCase.Mutations[ line[0] ].Raw.add( ( line[8] , line[36] ) )        
# ======================================================================================================

# ======================================================================================================
def AddRnaSeqFileToCase( aCase , aFile ):
  lStarCounts = StarCounts()
  
  for line in utf8reader( aFile ):
    if not line.startswith( "ENSG" ): continue
    line = [ i.strip() for i in line.split( "\t" , maxsplit = 7 ) ]
    lGeneName , lGeneType , lValue = line[1] , line[2] , float( line[6] ) # Gene-name , Gene-type , TPM_unstranded
    if True: #lGeneType == "protein_coding" : 
      if not lGeneName in StarCounts.GeneCatalogue:
        StarCounts.GeneCatalogue[ lGeneName ] = len( StarCounts.GeneCatalogue )
        lStarCounts.Genes.append( lValue )
      else:
        lStarCounts.Genes[ StarCounts.GeneCatalogue[lGeneName] ] = lValue
        
  aCase.StarCounts.append( lStarCounts )
# ======================================================================================================



# ======================================================================================================
print( "Creating Cases" , flush=True )
lInfo = [ "cases.case_id" , "cases.diagnoses.age_at_diagnosis" , "cases.project.disease_type" , "cases.project.primary_site" , "experimental_strategy" ]       
lFileInfo = GetJson( Or( And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "RNA-Seq" ), Eq( "files.analysis.workflow_type" , "STAR - Counts"), Eq( "files.data_type" , "Gene Expression Quantification"), Eq( "files.access" , "open") ) ,       
                         And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "WXS" ), Eq( "files.access" , "open") ) ) , lInfo )

lCases = {}

length = len( lFileInfo )
for i in tqdm.tqdm( range( 0 , length , 50 ) , ncols=Ncol , desc="Getting Data" ):  
  lSlice = lFileInfo[ i : min( length , i + 50 ) ]
  lFileIds = {}

  for j in lSlice:  
    CaseId = j["cases"][0]["case_id"]
    FileId = j["id"]
    ExperimentalStrategy = j["experimental_strategy"]

    lFileIds[ FileId ] = ( CaseId , ExperimentalStrategy )
    
    if not CaseId in lCases:
      DiseaseType = j["cases"][0]["project"]["disease_type"] 
      PrimarySite = j["cases"][0]["project"]["primary_site"]
      try :    AgeAtDiagnosis = int( j["cases"][0]["diagnoses"]["age_at_diagnosis"] )
      except : AgeAtDiagnosis = None   
      
      lCases[ CaseId ] = Case( CaseId , AgeAtDiagnosis , DiseaseType , PrimarySite )

  lTar = CachedGetTar( list( lFileIds.keys() ) )
    
  for lName in tqdm.tqdm( lTar.getmembers(), leave=False , ncols=Ncol ):                                    # Iterate over the files in the tarball
    if lName.name == "MANIFEST.txt": continue
    lFileId , _ = lName.name.split( "/" , maxsplit=1 )
    lCaseId , lStrategy  = lFileIds[ lFileId ]
    if lStrategy == "WXS" : AddWxsFileToCase(    lCases[ lCaseId ] , lTar.extractfile( lName ) )
    else:                   AddRnaSeqFileToCase( lCases[ lCaseId ] , lTar.extractfile( lName ) )

  lTar.close()

for i,j in tqdm.tqdm( lCases.items() , ncols=Ncol , desc="Classifying mutations" ):
  for k,l in j.Mutations.items():
    l.classify()
    
SaveCases( args.dest , lCases )
# ======================================================================================================