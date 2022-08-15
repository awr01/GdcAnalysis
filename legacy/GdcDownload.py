import requests, json, datetime, tarfile, io, codecs, gzip

# ======================================================================================================
# Magic to simplify writing filters
def And( *args ): return { "op":"and" , "content":[ i for i in args ] }
def Or( *args ):  return { "op":"or" , "content":[ i for i in args ] }
def Eq( key , value ):  return { "op":"=" , "content":{ "field" : key , "value": value } }
# ======================================================================================================

lFilter = Or(
            And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "RNA-Seq" ), Eq( "files.analysis.workflow_type" , "STAR - Counts"), Eq( "files.data_type" , "Gene Expression Quantification") ),
            And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "WXS" ), Eq( "files.access" , "open") )
         )

lTimeStamp = datetime.datetime.now().strftime( "%Y-%m-%d-%H-%M-%S" )
lDestName = f"GDC-{lTimeStamp}.tgz"
with tarfile.open( lDestName , mode='w|gz' ) as dest:

  # ----------------------------------------------------------------------------------------------------
  print( "Getting file lists" , flush=True )           
  lParams = { "filters":json.dumps( lFilter ) , "fields":"experimental_strategy,cases.case_id" , "size": "9999999" }                   
  lResponse = json.loads( requests.get( "https://api.gdc.cancer.gov/files" , params = lParams ).content.decode("utf-8") ) # Get response and convert to python
  lFileList = { i["id"] : ( i["cases"][0]["case_id"] , i["experimental_strategy"] ) for i in lResponse["data"]["hits"] }
  lFileIds = list( lFileList.keys() )
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( f"Downloading {len(lFileIds)} files as single epic tarball" , flush=True )
  lTarball = requests.post( "https://api.gdc.cancer.gov/data" , data = json.dumps( { "ids":lFileIds } ) , headers = { "Content-Type" : "application/json" } ).content
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( f"Reformatting GDC tarball to hierarchical format file '{lDestName}'" , flush=True )  
  with tarfile.open( fileobj = io.BytesIO( lTarball ) , mode='r' ) as src:
    lMembers = src.getmembers()
    for lName in lMembers:
      print( f"." , end="" , flush=True )    
      if lName.name == "MANIFEST.txt" : continue
      
      lFileId , lFileName = lName.name.split( "/" , maxsplit=1 )
      lCaseId , lMethod = lFileList[ lFileId ]
     
      lData = io.BytesIO( src.extractfile( lName ).read() )

      lInfo = tarfile.TarInfo( f'GDC-{lTimeStamp}/{lCaseId}/{lMethod}/{lFileName}' )
      lInfo.size = lData.getbuffer().nbytes
      dest.addfile( lInfo , lData )
  print( flush = True ) 
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( "Complete" , flush = True )  
  # ----------------------------------------------------------------------------------------------------

           
