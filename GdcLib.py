import tarfile, tqdm , io , _pickle , multiprocessing , os , gzip

Ncol = 150

# ======================================================================================================
class Meta(type):
  def __repr__(cls): return getattr( cls, 'class_str' )

# Mutation enumerations
class WildType(metaclass=Meta):       class_str = "wild-type" 
class SingleMutation(metaclass=Meta): class_str = "single-mutation"
class MultiMutation(metaclass=Meta):  class_str = "multi-mutation"
class SilentOrSplice(metaclass=Meta): class_str = "silent-or-splice"

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
class GeneCatalogueEntry:
  def __init__( self , index , typename ):
    self.index = index
    self.type = typename

class StarCounts:
  GeneCatalogue = {}

  def __init__( self ): 
    self.TpmUnstranded = [ None for i in range( len( StarCounts.GeneCatalogue ) ) ]

  def __getitem__( self , aGene ): 
    return self.TpmUnstranded[ StarCounts.GeneCatalogue[ aGene ].index ]
# ======================================================================================================

# ======================================================================================================
class Case:
  def __init__( self , CaseId , AgeAtDiagnosis , DiseaseType ):
    self.CaseId = CaseId
    self.AgeAtDiagnosis = AgeAtDiagnosis
    self.DiseaseType = DiseaseType

    self.Mutations = {}
    self.StarCounts = [] 
    
  def GetMutations( self , aGene , default = None ):
    if aGene in self.Mutations: return self.Mutations[ aGene ]
    return default
    
  def FormattedAge( self ):
    years = int( self.AgeAtDiagnosis / 365.25 )
    days = int( self.AgeAtDiagnosis - ( 365.25 * years ) )
    return f"{years:02}|{days:03}"    
# ======================================================================================================

# ======================================================================================================
def DumpToTar( dest , obj , name ):
  lData = io.BytesIO( _pickle.dumps( obj ) )    
  lInfo = tarfile.TarInfo( name )
  lInfo.size = lData.getbuffer().nbytes
  dest.addfile( lInfo , lData )

def SaveCases( aFilename , aCases ):
  with tarfile.open( aFilename , mode='a' ) as dest:
    DumpToTar( dest , StarCounts.GeneCatalogue , "@GeneCatalogue" )
    for i in tqdm.tqdm( aCases , ncols=Ncol , desc="Saving to disk" ): DumpToTar( dest , i , i.CaseId )
# ======================================================================================================

# # ======================================================================================================
# def LoadCases( aFilename ):
#   lCases = []
#   print( f"Opening '{aFilename}'" , flush=True )
#   with tarfile.open( aFilename , mode = 'r' ) as src:
#     StarCounts.GeneCatalogue = _pickle.loads( src.extractfile( "@GeneCatalogue" ).read() )       

#     for lName in tqdm.tqdm( src.getmembers() , ncols=Ncol , desc="Loading cases" ):
#       if lName.name[0] != "@": lCases.append( _pickle.loads( src.extractfile( lName ).read() ) )
      
#   return lCases
# # ======================================================================================================

# # ======================================================================================================
# def LoadCases( aFilename , aCaseIds ):
#   lCases = []
#   with tarfile.open( aFilename , mode = 'r' ) as src:
#     StarCounts.GeneCatalogue = _pickle.loads( src.extractfile( "@GeneCatalogue" ).read() )       

#     for lName in tqdm.tqdm( aCaseIds , leave=False , ncols=Ncol , desc="Loading cases" ):
#       lCases.append( _pickle.loads( src.extractfile( lName ).read() ) )
      
#   return lCases
# # ======================================================================================================

# # ======================================================================================================
# def LoadAndForEach( aFilename , aFn , Before = None , After = None ):
#   print( f"Opening '{aFilename}'" , flush=True )

#   with tarfile.open( aFilename , mode = 'r' ) as src:
#     StarCounts.GeneCatalogue = _pickle.loads( src.extractfile( "@GeneCatalogue" ).read() )       
#     if not Before is None:  Before()    
#     for lName in tqdm.tqdm( src.getmembers() , ncols=Ncol , desc="Load and analyze" ):
#       if lName.name[0] != "@" : aFn( _pickle.loads( src.extractfile( lName ).read() ) )      
#     if not After is None:  After()
# # ======================================================================================================


# ======================================================================================================
def ClassifyHandler( Args ):
  aFilename , ForEachClass , Class , Ids = Args
  index = multiprocessing.current_process()._identity[0]
  with tarfile.open( aFilename , mode = 'r' ) as src:    
    Cases = [ _pickle.loads( src.extractfile( lName ).read() ) for lName in tqdm.tqdm( Ids , leave=False , ncols=Ncol , desc=f"{Class}: Loading" , position=index ) ]
    return Class ,  ForEachClass( Class , Cases , index )

def LoadAndClassify( aFilename , ClassifyFn , ForEachClassFn , FinallyFn , maxthreads=None , cachefile=None ):

  # -------------------------
  # If a cache exists, use it
  if (not cachefile is None) and os.path.isfile( cachefile ):
    print( f"Loading cache: '{cachefile}': If you have recently downloaded new data, you may need to delete this file.")
    with gzip.open( cachefile , 'rb' ) as src: Data = _pickle.load( src )  
  # -------------------------
  # Else analyze the code
  else:
    with tarfile.open( aFilename , mode = 'r' ) as src:
      StarCounts.GeneCatalogue = _pickle.loads( src.extractfile( "@GeneCatalogue" ).read() )       

      Classes = {}
      for lName in tqdm.tqdm( src.getmembers() , leave=False , ncols=Ncol , desc=f"Load and classify {aFilename}" ):
        if lName.name[0] == "@" : continue
        Class = ClassifyFn( _pickle.loads( src.extractfile( lName ).read() ) )      
        if Class in Classes: Classes[ Class ].append( lName )
        else:                Classes[ Class ] = [ lName ] 
      
    with multiprocessing.Pool( maxthreads ) as pool:
      Data = {}
      generator = pool.imap( ClassifyHandler , (( aFilename , ForEachClassFn , k , v ) for k,v in sorted( Classes.items()) ) ) # Bodgy hack to bypass the limitations of imap
      for res in tqdm.tqdm( generator , ncols=Ncol, desc="Analysing", total=len(Classes) ): Data[ res[0] ] = res[1]

    if (not cachefile is None):
      print( f"Writing cache: '{cachefile}'" )
      with gzip.open( cachefile , 'wb' ) as dest: _pickle.dump( Data , dest )
  # -------------------------

  print( f"Finally" )
  FinallyFn( Data )
# ======================================================================================================