from GdcLib import *
from GdcUtils import *

# ======================================================================================================
def Init():
  dest.write( "CaseId\tAge at diagnosis\tDisease type" )                                                         # Write the header-line to the destination tsv
  for i in args.mutations: dest.write( "\tMutations-" + i )
  for i in args.genes:     dest.write( "\tStarCount-" + i )
  dest.write( "\n" )
# ======================================================================================================

# ======================================================================================================
def Analyze( aCase ):
  for lStarCount in aCase.StarCounts:                                                                            # Iterate over each star-count file
    dest.write( f"{aCase.CaseId}\t{aCase.FormattedAge()}\t{aCase.DiseaseType}" )                                 # Write the Case info
    for i in args.mutations: dest.write( "\t" + str( aCase.GetMutations( i , "" ) ) )                            # Write the mutations for the requested genes
    for i in args.genes:     dest.write( "\t" + str( lStarCount[i] ) )                                           # Write the star-counts for the requested genes
    dest.write( "\n" )                                                                                           # Write a new-line
# ======================================================================================================

# ======================================================================================================
with open( args.dest , "w" ) as dest: LoadAndForEach( args.src , Analyze , Init ) # Open the destination tsv and call the file-handle 'dest', then load src and analyze on-the-fly, making use of the optional Before argument
# ======================================================================================================
