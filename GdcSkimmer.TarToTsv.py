from GdcLib import *
from GdcUtils import *

# ======================================================================================================
def Init():
  dest.write( "CaseId" )                                                         # Write the header-line to the destination tsv
  for i in args.mutations: dest.write( "\tMutations-" + i )
  for i in args.genes:     dest.write( "\tStarCount-" + i )
  dest.write( "\n" )
# ======================================================================================================

# ======================================================================================================
def Analyze( aCase ):
  for lStarCount in aCase.StarCounts:                                            # Then iterate over each star-count file
    dest.write( aCase.CaseId )
    for i in args.mutations: dest.write( "\t" + str( aCase.GetMutations( i , "" ) ) )         
    for i in args.genes:     dest.write( "\t" + str( lStarCount[i] ) )
    dest.write( "\n" ) 
# ======================================================================================================

# ======================================================================================================
with open( args.dest , "w" ) as dest: LoadAndForEach( args.src , Analyze , Init ) # Open the destination tsv and call the file-handle 'dest', then load src and analyze on-the-fly, making use of the optional Before argument
# ======================================================================================================
