import requests, json, datetime, tarfile, io, codecs, gzip

TimeStamp = datetime.datetime.now().strftime( "%Y-%m-%d-%H-%M-%S" )
 
# ==============================================================================================================================================
def GetFileList( aFilter ):
  print( "Get file lists" , flush=True )           
  lParams = { "filters":json.dumps( aFilter ) , "fields":"experimental_strategy,cases.case_id" , "size": "9999999" }                   
  lResponse = json.loads( requests.get( "https://api.gdc.cancer.gov/files" , params = lParams ).content.decode("utf-8") ) # Get response and convert to python
  return { i["id"] : ( i["cases"][0]["case_id"] , i["experimental_strategy"] ) for i in lResponse["data"]["hits"] }
# ==============================================================================================================================================  

# ==============================================================================================================================================
def GetFiles( aFileList , aRawFileName = None ):
  print( f"Downloading {len(aFileList)} files as single epic tarball" , flush=True )
  lTarball = requests.post( "https://api.gdc.cancer.gov/data" , data = json.dumps( { "ids":aFileList } ) , headers = { "Content-Type" : "application/json" } ).content

  if not aRawFileName is None:
    print( f"Save the raw tarball to file '{aRawFileName}'" , flush=True )
    with open( aRawFileName , "wb") as lOutFile: lOutFile.write( aData )
  
  return lTarball
# ==============================================================================================================================================
  
# ==============================================================================================================================================
def SaveReformattedTarball( aFileList , aTarball , aFileName = f"GDC-{TimeStamp}.tgz" ):  
  print( f"Reformatting GDC tarball to hierarchical format file '{aFileName}'" , flush=True )
  
  with tarfile.open( fileobj = io.BytesIO( aTarball ) , mode='r' ) as src , tarfile.open( aFileName , mode='w|gz' ) as dest:
    lMembers = src.getmembers()
    count , total = 1 , len(lMembers)
    for lName in lMembers:
      print( f"{count:4}/{total:4}\r" , end="" , flush=True )
      count += 1      
      if lName.name == "MANIFEST.txt" : continue
      
      lFileId , lFileName = lName.name.split( "/" , maxsplit=1 )
      lCaseId , lMethod = aFileList[ lFileId ]
     
      lData = io.BytesIO( src.extractfile( lName ).read() )

      lInfo = tarfile.TarInfo( f'GDC-{TimeStamp}/{lCaseId}/{lMethod}/{lFileName}' )
      lInfo.size = lData.getbuffer().nbytes
      dest.addfile( lInfo , lData )
  
  return aFileName
# ==============================================================================================================================================

# ==============================================================================================================================================
def And( *args ): return { "op":"and" , "content":[ i for i in args ] }
def Or( *args ):  return { "op":"or" , "content":[ i for i in args ] }
def Eq( key , value ):  return { "op":"=" , "content":{ "field" : key , "value": value } }
# ==============================================================================================================================================


# ==============================================================================================================================================
class Case:

  class Mutation:
    def __init__( self , Type , Details ): self.Type , self.Details = Type , Details
    def __str__( self ): return ( f"{self.Type}[{self.Details}]" if len( self.Details ) else self.Type )

  def __init__( self , aCaseId ):
    self.CaseId = aCaseId
    self.Mutations = None
    self.StarCounts = None
  
  def getMutations( self , aGene ):
    if self.Mutations is None : return None
    return [ i[1] for i in self.Mutations if i[0] == aGene ]
# ==============================================================================================================================================


# ==============================================================================================================================================
def LoadFormattedTarball( aFileName , aGeneList = None ):
  print( f"Opening formatted tarball '{aFileName}'" , flush=True )
  with tarfile.open( aFileName ) as src:

    utf8reader = codecs.getreader( 'utf-8' )

    lCases = {}
    lMembers = src.getmembers()
    count , total = 1 , len(lMembers)
    
    for lName in lMembers: 
      print( f"{count:4}/{total:4}\r" , end="" , flush=True )
      count += 1
      # if count == 10: break
      _ , lCaseId , lMethod , _ = lName.name.split( "/" )     
      
      if not lCaseId in lCases : lCases[ lCaseId ] = Case( lCaseId )
      lCase = lCases[ lCaseId ]

      # ------------------------------------------------------------------------------
      if lMethod == 'WXS': 
        if lCase.Mutations is None: lCase.Mutations = []    
        for line in utf8reader( gzip.GzipFile( fileobj = src.extractfile( lName ) ) ):
          if line[0] == "#" or line.startswith( "Hugo_Symbol" ) : continue
          line = [ i.strip() for i in line.split( "\t" ) ]
          Name , Type , Details = line[0] , line[8] , line[36]
          lCase.Mutations.append( ( Name , Case.Mutation( Type , Details ) ) )     
      # ------------------------------------------------------------------------------      
      else:                
        if lCase.StarCounts is None: lCase.StarCounts = []  
        lData = {}
        for line in utf8reader( src.extractfile( lName ) ):
          line = [ i.strip() for i in line.split( "\t" ) ]
          if not line[0].startswith( "ENSG" ): continue
          lGeneName , lValue = line[1] , float( line[6] ) # Gene , TPM_unstranded
          if ( aGeneList is None ) or ( lGeneName in aGeneList ) : lData[ lGeneName ] = lValue
        lCase.StarCounts.append( lData )      
      # ------------------------------------------------------------------------------
      
    print( "=== Complete ===" , flush=True )  
    return lCases
# ==============================================================================================================================================

