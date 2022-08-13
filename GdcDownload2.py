import argparse , requests, json, tarfile, io, codecs, gzip , tqdm , os , hashlib , bz2 , _pickle

# Read the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--prefix' , required=True , help='The prefix that will be used for all files')
args = parser.parse_args()

utf8reader = codecs.getreader( 'utf-8' )


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
  lFile = f".cache.{lParams}.pkl.bz"
  
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
  lFile = f".cache.{lKeys}.tgz"

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
class tWildType:
  def __repr__( self ) : return "wild-type"
  
class tSingleMutation:
  def __repr__( self ) : return "single-mutation"
  
class tMultiMutation:
  def __repr__( self ) : return "multi-mutation"
  
class tSilentOrSplice:
  def __repr__( self ) : return "silent-or-splice"

WildType = tWildType()
SingleMutation = tSingleMutation()
MultiMutation = tMultiMutation()
SilentOrSplice = tSilentOrSplice()

class WxsFile:

  class Mutation:  
    def __init__( self ):
      self.Raw = set()
      self.Classification = None
      self.IncludesSilents = False
      self.IncludesSplices = False
      
    def classify( self ):
      lMut = list( self.Raw )
      
      l0 = len( lMut )
      if l0 == 0:                             # If there are no mutations
        self.Classification = WildType
        return

      lMut = list( filter( lambda a: ("Splice" not in a[0]) , lMut ) ) # Remove all splices from the mutation list
      l1 = len( lMut )      
      if l1 != l0: self.IncludesSplices = True

      lMut = list( filter( lambda a: (a[0] != "Silent") , lMut ) )     # Remove all silent from the mutation list
      l2 = len( lMut )      
      if l2 != l1: self.IncludesSilents = True

      if   l2 == 0: self.Classification = SilentOrSplice       
      elif l2 == 1: self.Classification = SingleMutation
      else:         self.Classification = MultiMutation
      
    def __str__( self ):
      return f"{str(self.Classification):20}\t{self.IncludesSilents}\t{self.IncludesSplices}\t{self.Raw}"

  def __init__( self , lWxsFile ):  
    self.Mutations = {}
  
    for line in utf8reader( gzip.GzipFile( fileobj = lWxsFile ) ): # Iterate over each line in the file
      if line[0] == "#" or line.startswith( "Hugo_Symbol" ) : continue             # Ignore comments and headers
      line = [ i.strip() for i in line.split( "\t" , maxsplit = 37 ) ]              # Split the line at tabs up to where we need it      
      if not line[0] in self.Mutations: self.Mutations[ line[0] ] = WxsFile.Mutation()
      self.Mutations[ line[0] ].Raw.add( ( line[8] , line[36] ) )        
      
    for i,j in self.Mutations.items(): j.classify()
# ======================================================================================================


# ======================================================================================================
class RnaSeqFile:
  def __init__( self , lRnaSeqFile ): 
    self.Genes = {}

    for line in utf8reader( lRnaSeqFile ):
      if not line.startswith( "ENSG" ): continue
      line = [ i.strip() for i in line.split( "\t" , maxsplit = 7 ) ]
      lGeneName , lGeneType , lValue = line[1] , line[2] , float( line[6] ) # Gene-name , Gene-type , TPM_unstranded
      if lGeneType == "protein_coding" : self.Genes[ lGeneName ] = lValue
# ======================================================================================================


# ======================================================================================================
class Case:
  def __init__( self , CaseId , data ):
    self.CaseId = CaseId
    self.WxsFileIds = { data["id"] : None }
    self.RnaSeqFileIds = {}
    
    try :    self.AgeAtDiagnosis = int( data["cases"][0]["diagnoses"]["age_at_diagnosis"] )
    except : self.AgeAtDiagnosis = None  
    try:     self.DiseaseType = data["cases"][0]["project"]["disease_type"] 
    except : self.DiseaseType = None    
    try:     self.PrimarySite = data["cases"][0]["project"]["primary_site"]
    except:  self.PrimarySite = None
# ======================================================================================================









# ======================================================================================================
print( "Creating Catalogue" , flush=True )
lInfo    = [ "cases.case_id" ]
lAuxInfo = [ "cases.case_id" , "cases.diagnoses.age_at_diagnosis" , "cases.project.disease_type" , "cases.project.primary_site" ]       
lRnaSeq  = GetJson( And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "RNA-Seq" ), Eq( "files.analysis.workflow_type" , "STAR - Counts"), Eq( "files.data_type" , "Gene Expression Quantification"), Eq( "files.access" , "open") ) , lInfo )        
lWxs     = GetJson( And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "WXS" ), Eq( "files.access" , "open") ) , lAuxInfo )
 
lCatalogue = {}
for i in lWxs:
  CaseId = i["cases"][0]["case_id"]
  if not CaseId in lCatalogue : 
    lCatalogue[ CaseId ] = Case( CaseId , i )
  else:                         
    lCatalogue[ CaseId ].WxsFileIds[ i["id"] ] = None

for i in lRnaSeq:
  CaseId = i["cases"][0]["case_id"]
  if not CaseId in lCatalogue : 
    continue # We don't care about cases with no WXS
  else:
    lCatalogue[ CaseId ].RnaSeqFileIds[ i["id"] ] = None
  
for i in list( lCatalogue.keys() ):
  if len( lCatalogue[i].RnaSeqFileIds ) == 0 : del lCatalogue[i] # If there is no RnaSeq - delete that too 
# ======================================================================================================

# ======================================================================================================
print( "Getting WXS info" , flush=True )
lWxsFileIds = {}
for i,j in lCatalogue.items(): 
  for k in j.WxsFileIds: lWxsFileIds[ k ] = i
    
lWxsFiles = GetTar( list( lWxsFileIds.keys() ) , 100 )
# ======================================================================================================

# ======================================================================================================
print( "Extracting mutations" , flush=True )
for lName in tqdm.tqdm( lWxsFiles.getmembers() , ncols=100 ):                                    # Iterate over the files in the tarball
  lFileId , lFileName = lName.name.split( "/" , maxsplit=1 )
  lCaseId = lWxsFileIds[ lFileId ]
  lCatalogue[ lCaseId ].WxsFileIds[ lFileId ] = WxsFile( lWxsFiles.extractfile( lName ) )
# ======================================================================================================
  
# ======================================================================================================
print( "Getting RnaSeq info" , flush=True )
lRnaSeqFileIds = {}
for i,j in lCatalogue.items(): 
  for k in j.RnaSeqFileIds: lRnaSeqFileIds[ k ] = i
    
lRnaSeqFiles = GetTar( list( lRnaSeqFileIds.keys() ) , 100 )
# ======================================================================================================    

# ======================================================================================================
print( "Extracting tpm-unstranded" , flush=True )
for lName in tqdm.tqdm( lRnaSeqFiles.getmembers() , ncols=100 ):                                    # Iterate over the files in the tarball
  lFileId , lFileName = lName.name.split( "/" , maxsplit=1 )
  lCaseId = lRnaSeqFileIds[ lFileId ]
  lCatalogue[ lCaseId ].RnaSeqFileIds[ lFileId ] = RnaSeqFile( lRnaSeqFiles.extractfile( lName ) )
# ======================================================================================================

# ======================================================================================================
with bz2.BZ2File( f"{args.prefix}.pkl.bz2" , 'w' ) as dest: _pickle.dump( lCatalogue , dest )
# ======================================================================================================
