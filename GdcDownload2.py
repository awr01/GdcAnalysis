import argparse , requests, json, tarfile, io, codecs, gzip , tqdm , os , hashlib , bz2 , _pickle

from GdcLib2 import *

# Read the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--prefix' , required=True , help='The prefix that will be used for all files')
args = parser.parse_args()


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
def GetTar( aKeys , aBatchsize ):
  lKeys = hashlib.md5( str( sorted( aKeys ) ).encode('utf-8') ).hexdigest()
  lFile = f".cache/{lKeys}.tgz"

  if os.path.isfile( lFile ): 
    print( "> Using cached tar" , flush=True )
  else: 
    print( "> Downloading tar" , flush=True )
    with tarfile.open( lFile , mode='w|gz' ) as dest:
      for i in tqdm.tqdm( range( 0 , len( aKeys ) , aBatchsize ) , ncols=100 ):
        lSubsetKeys = aKeys[ i : min( i+aBatchsize , len( aKeys ) ) ]

        lResponse = requests.post( "https://api.gdc.cancer.gov/data" , data = json.dumps( { "ids":lSubsetKeys } ) , headers = { "Content-Type" : "application/json" } )
        if lResponse.status_code != 200: raise Exception( f"Requests returned status-code {lResponse.status_code}" )  
       
        with tarfile.open( fileobj = io.BytesIO( lResponse.content ) , mode='r' ) as src:
          for lName in src.getmembers():  
            if lName.name == "MANIFEST.txt" : continue           
            lData = io.BytesIO( src.extractfile( lName ).read() )
            lInfo = tarfile.TarInfo( lName.name )
            lInfo.size = lData.getbuffer().nbytes
            dest.addfile( lInfo , lData )
          
  src = tarfile.open( lFile , mode='r:gz' )
  return src
# ======================================================================================================










# ======================================================================================================
print( "Creating Cases" , flush=True )
lInfo    = [ "cases.case_id" ]
lAuxInfo = [ "cases.case_id" , "cases.diagnoses.age_at_diagnosis" , "cases.project.disease_type" , "cases.project.primary_site" ]       
lRnaSeq  = GetJson( And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "RNA-Seq" ), Eq( "files.analysis.workflow_type" , "STAR - Counts"), Eq( "files.data_type" , "Gene Expression Quantification"), Eq( "files.access" , "open") ) , lInfo )        
lWxs     = GetJson( And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "WXS" ), Eq( "files.access" , "open") ) , lAuxInfo )
 
lCases = {}
for i in lWxs:
  CaseId = i["cases"][0]["case_id"]
  if not CaseId in lCases : 
    lCases[ CaseId ] = Case( CaseId , i )
  else:                         
    lCases[ CaseId ].WxsFileIds[ i["id"] ] = None

for i in lRnaSeq:
  CaseId = i["cases"][0]["case_id"]
  if not CaseId in lCases : 
    continue # We don't care about cases with no WXS
  else:
    lCases[ CaseId ].RnaSeqFileIds[ i["id"] ] = None
  
for i in list( lCases.keys() ):
  if len( lCases[i].RnaSeqFileIds ) == 0 : del lCases[i] # If there is no RnaSeq - delete that too 
# ======================================================================================================

# ======================================================================================================
print( "Getting WXS info" , flush=True )
lWxsFileIds = {}
for i,j in lCases.items(): 
  for k in j.WxsFileIds: lWxsFileIds[ k ] = i
    
lWxsFiles = GetTar( list( lWxsFileIds.keys() ) , 100 )
# ======================================================================================================

# ======================================================================================================
print( "Extracting mutations" , flush=True )
for lName in tqdm.tqdm( lWxsFiles.getmembers() , ncols=100 ):                                    # Iterate over the files in the tarball
  lFileId , lFileName = lName.name.split( "/" , maxsplit=1 )
  lCaseId = lWxsFileIds[ lFileId ]
  lCases[ lCaseId ].WxsFileIds[ lFileId ] = WxsFile( lWxsFiles.extractfile( lName ) )
# ======================================================================================================
  
# ======================================================================================================
print( "Getting RnaSeq info" , flush=True )
lRnaSeqFileIds = {}
for i,j in lCases.items(): 
  for k in j.RnaSeqFileIds: lRnaSeqFileIds[ k ] = i
    
lRnaSeqFiles = GetTar( list( lRnaSeqFileIds.keys() ) , 100 )
# ======================================================================================================    

# ======================================================================================================
print( "Extracting tpm-unstranded" , flush=True )
for lName in tqdm.tqdm( lRnaSeqFiles.getmembers() , ncols=100 ):                                    # Iterate over the files in the tarball
  lFileId , lFileName = lName.name.split( "/" , maxsplit=1 )
  lCaseId = lRnaSeqFileIds[ lFileId ]
  lCases[ lCaseId ].RnaSeqFileIds[ lFileId ] = RnaSeqFile( lRnaSeqFiles.extractfile( lName ) )
# ======================================================================================================


# ======================================================================================================
SaveCases( f"{args.prefix}.tgz" , lCases )
# ======================================================================================================
