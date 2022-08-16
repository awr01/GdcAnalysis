import argparse , tarfile , codecs , gzip , tqdm
from numpy import mean , var , warnings
warnings.filterwarnings('ignore')


def Classify( Mutations ):
  if len( Mutations ) == 0: return 0                                      # If there are no mutations

  Mutations = list( filter( lambda a: ("Splice" not in a) , Mutations ) ) # Remove all splices from the mutation list
  if len( Mutations ) == 0: return 2                                      # If there are now no mutations (i.e. there were only slices)

  Mutations = list( filter( lambda a: (a != "Silent") , Mutations ) )     # Remove all silent from the mutation list
  if len( Mutations ) == 0: return 3                                      # If there are now no mutations (i.e. there were only silents)           

  return 1 # There were other mutations
 



utf8reader = codecs.getreader( 'utf-8' )

# Read the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')
args = parser.parse_args()


Summary = {}  # Initialize empty lists to store case-ids for the different mutation-classifications

print( f"Opening formatted tarball '{args.src}'" , flush=True )
with tarfile.open( args.src ) as src , open( args.dest , "w" ) as dest: # Open the formatted tarball file and the destination files
 
  # ----------------------------------------------------------------------------------------------------  
  print( "Getting mutations" , flush = True )
  for lName in tqdm.tqdm( src.getmembers() , ncols=100 ):                                    # Iterate over the files in the tarball
    _ , lCaseId , lMethod , _ = lName.name.split( "/" )     # Extract the case-id and file-type from the name     
    if lMethod == "WXS":                                    # Only look at WXS files
      Mutations = { "ATRX" : set() }                                        # For each file initialize an empty mutations list
      
      for line in utf8reader( gzip.GzipFile( fileobj = src.extractfile( lName ) ) ): # Iterate over each line in the file
        if line[0] == "#" or line.startswith( "Hugo_Symbol" ) : continue             # Ignore comments and headers
        line = [ i.strip() for i in line.split( "\t" , maxsplit = 9 ) ]              # Split the line at tabs up to where we need it       
        if not line[0] in Mutations : Mutations[ line[0] ] = set()        
        Mutations[ line[0] ].add( line[8] )                                          # If the mutation is on ATRX, add it to the list

      ATRX = Classify( Mutations[ "ATRX" ] )

      for i,j in Mutations.items():
        if i == "ATRX" : continue        
        if not i in Summary: Summary[ i ] = [ [ 0 for y in range(4)] for x in range(4) ]
        Summary[ i ][ ATRX ][ Classify( j ) ] += 1

  dest.write( "\tATRX-WT\t\t\t\tATRX-Mut\t\t\t\tATRX-splice\t\t\t\tATRX-silent\t\t\t\t\n" )
  dest.write( "Gene\tGene-WT\tGene-Mut\tGene-splice\tGene-silent\tGene-WT\tGene-Mut\tGene-splice\tGene-silent\tGene-WT\tGene-Mut\tGene-splice\tGene-silent\tGene-WT\tGene-Mut\tGene-splice\tGene-silent\n" )
  for i,j in Summary.items():
    dest.write( i )
    for x in range( 4 ):
      for y in range( 4 ):
        dest.write( f"\t{j[x][y]}" )
    dest.write( "\n" )

  # print( flush = True ) 
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( "Complete" , flush = True )  
  # ----------------------------------------------------------------------------------------------------
