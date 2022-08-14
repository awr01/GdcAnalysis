import argparse , requests, json, tarfile, io, codecs, gzip , tqdm , os , hashlib , bz2 , _pickle

utf8reader = codecs.getreader( 'utf-8' )

GeneCatalogue = []


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
    self.Genes = [ None for i in range( len( GeneCatalogue ) ) ]

    cnt = 0
    for line in utf8reader( lRnaSeqFile ):
      if not line.startswith( "ENSG" ): continue
      line = [ i.strip() for i in line.split( "\t" , maxsplit = 7 ) ]
      lGeneName , lGeneType , lValue = line[1] , line[2] , float( line[6] ) # Gene-name , Gene-type , TPM_unstranded
      if lGeneType == "protein_coding" : 
        try:
          if GeneCatalogue[ cnt ] == lGeneName:
            self.Genes[ cnt ] = lValue
          else:
            lGeneIndex = GeneCatalogue.index( lGeneName )
            self.Genes[ lGeneIndex ] = lValue
        except:
          GeneCatalogue.append( lGeneName )
          self.Genes.append( lValue )
        cnt += 1
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
def SaveCases( aFilename , aCases ):
  with tarfile.open( aFilename , mode='w|gz' ) as dest:

    lData = io.BytesIO( _pickle.dumps( GeneCatalogue ) )    
    lInfo = tarfile.TarInfo( "GeneCatalogue" )
    lInfo.size = lData.getbuffer().nbytes
    dest.addfile( lInfo , lData )
    
    for i,j in aCases.items():
      lData = io.BytesIO( _pickle.dumps( j ) )    
      lInfo = tarfile.TarInfo( i )
      lInfo.size = lData.getbuffer().nbytes
      dest.addfile( lInfo , lData )
# ======================================================================================================

# ======================================================================================================
def LoadCases( aFilename ):
  print( "Loading cases" )
  lCases = {}
  with tarfile.open( aFilename , mode='r:gz' ) as src:
    for lName in tqdm.tqdm( src.getmembers() , ncols=100 ):
      if lName.name == "GeneCatalogue" : GeneCatalogue        = _pickle.loads( src.extractfile( lName ).read() )
      else:                              lCases[ lName.name ] = _pickle.loads( src.extractfile( lName ).read() )
  return lCases
# ======================================================================================================
