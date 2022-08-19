import tarfile, tqdm , io , _pickle

Ncol = 150

# ======================================================================================================
class Meta(type):
  def __repr__(cls): return getattr( cls, 'class_str' )

# Mutation enumerations
class WildType(metaclass=Meta):       class_str = "wild-type" 
class SingleMutation(metaclass=Meta): class_str = "single-mutation"
class MultiMutation(metaclass=Meta):  class_str = "multi-mutation"
class SilentOrSplice(metaclass=Meta): class_str = "silent-or-splice"

# Gene-type enumerations
class IG_C_gene(metaclass=Meta): class_str = "IG_C_gene" 
class IG_C_pseudogene(metaclass=Meta): class_str = "IG_C_pseudogene" 
class IG_D_gene(metaclass=Meta): class_str = "IG_D_gene" 
class IG_J_gene(metaclass=Meta): class_str = "IG_J_gene" 
class IG_J_pseudogene(metaclass=Meta): class_str = "IG_J_pseudogene" 
class IG_pseudogene(metaclass=Meta): class_str = "IG_pseudogene" 
class IG_V_gene(metaclass=Meta): class_str = "IG_V_gene" 
class IG_V_pseudogene(metaclass=Meta): class_str = "IG_V_pseudogene" 
class lncRNA(metaclass=Meta): class_str = "lncRNA" 
class miRNA(metaclass=Meta): class_str = "miRNA" 
class misc_RNA(metaclass=Meta): class_str = "misc_RNA" 
class Mt_rRNA(metaclass=Meta): class_str = "Mt_rRNA" 
class Mt_tRNA(metaclass=Meta): class_str = "Mt_tRNA" 
class polymorphic_pseudogene(metaclass=Meta): class_str = "polymorphic_pseudogene" 
class processed_pseudogene(metaclass=Meta): class_str = "processed_pseudogene" 
class protein_coding(metaclass=Meta): class_str = "protein_coding" 
class pseudogene(metaclass=Meta): class_str = "pseudogene" 
class ribozyme(metaclass=Meta): class_str = "ribozyme" 
class rRNA(metaclass=Meta): class_str = "rRNA" 
class rRNA_pseudogene(metaclass=Meta): class_str = "rRNA_pseudogene" 
class scaRNA(metaclass=Meta): class_str = "scaRNA" 
class scRNA(metaclass=Meta): class_str = "scRNA" 
class snoRNA(metaclass=Meta): class_str = "snoRNA" 
class snRNA(metaclass=Meta): class_str = "snRNA" 
class sRNA(metaclass=Meta): class_str = "sRNA" 
class TEC(metaclass=Meta): class_str = "TEC" 
class TR_C_gene(metaclass=Meta): class_str = "TR_C_gene" 
class TR_D_gene(metaclass=Meta): class_str = "TR_D_gene" 
class TR_J_gene(metaclass=Meta): class_str = "TR_J_gene" 
class TR_J_pseudogene(metaclass=Meta): class_str = "TR_J_pseudogene" 
class TR_V_gene(metaclass=Meta): class_str = "TR_V_gene" 
class TR_V_pseudogene(metaclass=Meta): class_str = "TR_V_pseudogene" 
class transcribed_processed_pseudogene(metaclass=Meta): class_str = "transcribed_processed_pseudogene" 
class transcribed_unitary_pseudogene(metaclass=Meta): class_str = "transcribed_unitary_pseudogene" 
class transcribed_unprocessed_pseudogene(metaclass=Meta): class_str = "transcribed_unprocessed_pseudogene" 
class translated_processed_pseudogene(metaclass=Meta): class_str = "translated_processed_pseudogene" 
class translated_unprocessed_pseudogene(metaclass=Meta): class_str = "translated_unprocessed_pseudogene" 
class unitary_pseudogene(metaclass=Meta): class_str = "unitary_pseudogene" 
class unprocessed_pseudogene(metaclass=Meta): class_str = "unprocessed_pseudogene" 
class vault_RNA(metaclass=Meta): class_str = "vault_RNA" 

  # # GeneTypeLUT = { "IG_C_gene":IG_C_gene , "IG_C_pseudogene":IG_C_pseudogene , "IG_D_gene":IG_D_gene , "IG_J_gene":IG_J_gene , "IG_J_pseudogene":IG_J_pseudogene , "IG_pseudogene":IG_pseudogene , 
                  # # "IG_V_gene":IG_V_gene , "IG_V_pseudogene":IG_V_pseudogene , "lncRNA":lncRNA , "miRNA":miRNA , "misc_RNA":misc_RNA , "Mt_rRNA":Mt_rRNA , "Mt_tRNA":Mt_tRNA , 
                  # # "polymorphic_pseudogene":polymorphic_pseudogene , "processed_pseudogene":processed_pseudogene , "protein_coding":protein_coding , "pseudogene":pseudogene , "ribozyme":ribozyme , 
                  # # "rRNA":rRNA , "rRNA_pseudogene":rRNA_pseudogene , "scaRNA":scaRNA , "scRNA":scRNA , "snoRNA":snoRNA , "snRNA":snRNA , "sRNA":sRNA , "TEC":TEC , "TR_C_gene":TR_C_gene , 
                  # # "TR_D_gene":TR_D_gene , "TR_J_gene":TR_J_gene , "TR_J_pseudogene":TR_J_pseudogene , "TR_V_gene":TR_V_gene , "TR_V_pseudogene":TR_V_pseudogene , 
                  # # "transcribed_processed_pseudogene":transcribed_processed_pseudogene , "transcribed_unitary_pseudogene":transcribed_unitary_pseudogene , 
                  # # "transcribed_unprocessed_pseudogene":transcribed_unprocessed_pseudogene , "translated_processed_pseudogene":translated_processed_pseudogene , 
                  # # "translated_unprocessed_pseudogene":translated_unprocessed_pseudogene , "unitary_pseudogene":unitary_pseudogene , "unprocessed_pseudogene":unprocessed_pseudogene , 
                  # # "vault_RNA":vault_RNA }
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
    
# ======================================================================================================



# ======================================================================================================
def DumpToTar( dest , obj , name ):
  lData = io.BytesIO( _pickle.dumps( obj ) )    
  lInfo = tarfile.TarInfo( name )
  lInfo.size = lData.getbuffer().nbytes
  dest.addfile( lInfo , lData )

def SaveCases( aFilename , aCases ):
  with tarfile.open( aFilename , mode='w|gz' ) as dest:
    DumpToTar( dest , StarCounts.GeneCatalogue , "@GeneCatalogue" )
    for i,j in tqdm.tqdm( aCases.items() , ncols=Ncol , desc="Saving to disk" ): DumpToTar( dest , j , i )
# ======================================================================================================

# ======================================================================================================
def LoadCases( aFilename ):
  lCases = {}
  with tarfile.open( aFilename , mode='r:gz' ) as src:
    StarCounts.GeneCatalogue = _pickle.loads( src.extractfile( "@GeneCatalogue" ).read() )       

    for lName in tqdm.tqdm( src.getmembers() , ncols=Ncol , desc="Loading cases" ):
      if lName.name[0] != "@": lCases[ lName.name ] = _pickle.loads( src.extractfile( lName ).read() )
      
  return lCases
# ======================================================================================================

# ======================================================================================================
def LoadAndForEach( aFilename , aFn , Before = None , After = None ):
  with tarfile.open( aFilename , mode='r:gz' ) as src:
    StarCounts.GeneCatalogue = _pickle.loads( src.extractfile( "@GeneCatalogue" ).read() )       
    if not Before is None:  Before()    
    for lName in tqdm.tqdm( src.getmembers() , ncols=Ncol , desc="Load and analyze" ):
      if lName.name[0] != "@" : aFn( _pickle.loads( src.extractfile( lName ).read() ) )      
    if not After is None:  After()
# ======================================================================================================
