import tarfile, tqdm , io , _pickle

Ncol = 150

# ======================================================================================================
class Meta(type):
  def __repr__(cls):
    return getattr( cls, 'class_str' )

class WildType(metaclass=Meta): class_str = "wild-type" 
class SingleMutation(metaclass=Meta): class_str = "single-mutation"
class MultiMutation(metaclass=Meta): class_str = "multi-mutation"
class SilentOrSplice(metaclass=Meta): class_str = "silent-or-splice"
# ======================================================================================================

# ======================================================================================================
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
    return " , ".join( [ f"{j[0]}({j[1]})" for j in self.Raw ] )
# ======================================================================================================

# ======================================================================================================
class StarCounts:
  GeneCatalogue = {}

  def __init__( self ): 
    self.Genes = [ None for i in range( len( StarCounts.GeneCatalogue ) ) ]

  def __getitem__( self , aGene ): 
    return self.Genes[ StarCounts.GeneCatalogue[ aGene ] ]
# ======================================================================================================

# ======================================================================================================
class Case:
  def __init__( self , CaseId , AgeAtDiagnosis , DiseaseType , PrimarySite ):
    self.CaseId = CaseId
    self.AgeAtDiagnosis = AgeAtDiagnosis
    self.DiseaseType = DiseaseType
    self.PrimarySite = PrimarySite

    self.Mutations = {}
    self.StarCounts = [] 
    
  def GetMutations( self , aGene , default = None ):
    if aGene in self.Mutations: return self.Mutations[ aGene ]
    return default
    
# ======================================================================================================



# ======================================================================================================
def SaveCases( aFilename , aCases ):
  with tarfile.open( aFilename , mode='w|gz' ) as dest:

    lData = io.BytesIO( _pickle.dumps( StarCounts.GeneCatalogue ) )    
    lInfo = tarfile.TarInfo( "GeneCatalogue" )
    lInfo.size = lData.getbuffer().nbytes
    dest.addfile( lInfo , lData )
    
    for i,j in tqdm.tqdm( aCases.items() , ncols=Ncol , desc="Saving to disk" ):
      lData = io.BytesIO( _pickle.dumps( j ) )    
      lInfo = tarfile.TarInfo( i )
      lInfo.size = lData.getbuffer().nbytes
      dest.addfile( lInfo , lData )
# ======================================================================================================

# ======================================================================================================
def LoadCases( aFilename ):
  lCases = {}
  with tarfile.open( aFilename , mode='r:gz' ) as src:
    for lName in tqdm.tqdm( src.getmembers() , ncols=Ncol , desc="Loading cases" ):
      if lName.name == "GeneCatalogue" : StarCounts.GeneCatalogue = _pickle.loads( src.extractfile( lName ).read() )
      else:                              lCases[ lName.name ]     = _pickle.loads( src.extractfile( lName ).read() )
      
  return lCases
# ======================================================================================================

# ======================================================================================================
def LoadAndForEach( aFilename , aFn , Before = None , After = None ):
  with tarfile.open( aFilename , mode='r:gz' ) as src:
    StarCounts.GeneCatalogue = _pickle.loads( src.extractfile( "GeneCatalogue" ).read() )       
    if not Before is None:  Before()    
    for lName in tqdm.tqdm( src.getmembers() , ncols=Ncol , desc="Load and analyze" ):
      if lName.name != "GeneCatalogue" : aFn( _pickle.loads( src.extractfile( lName ).read() ) )      
    if not After is None:  After()
# ======================================================================================================
