import argparse , tarfile , codecs , gzip
from scipy.stats import ttest_ind

parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')
args = parser.parse_args()

print( f"Opening formatted tarball '{args.src}'" , flush=True )
with tarfile.open( args.src ) as src , open( args.dest , "w" ) as dest , open( "raw-" + args.dest , "w" ) as dest2:

  utf8reader = codecs.getreader( 'utf-8' )
  lMembers = src.getmembers()
  
  # ----------------------------------------------------------------------------------------------------
  print( "Getting mutations" , flush = True )
  ATRXmut , ATRXwt , ATRXother = [] , [] , []
  for lName in lMembers: 
    _ , lCaseId , lMethod , _ = lName.name.split( "/" )  
     
    if lMethod == "WXS":
      Mutations = []
      
      for line in utf8reader( gzip.GzipFile( fileobj = src.extractfile( lName ) ) ):
        if line[0] == "#" or line.startswith( "Hugo_Symbol" ) : continue
        line = [ i.strip() for i in line.split( "\t" , maxsplit = 9 ) ]
        if line[0] == "ATRX" : Mutations.append( line[8] )          

      if len( Mutations ) == 0:
        ATRXwt.append( lCaseId )
        continue
      
      Mutations = list( filter( lambda a: (a != "Silent") and (a!="Splice_Site") , Mutations ) )
      if len( Mutations ) == 0:
        ATRXother.append( lCaseId )
        continue

      ATRXmut.append( lCaseId )
  # ----------------------------------------------------------------------------------------------------
     
  # ----------------------------------------------------------------------------------------------------
  print( "Getting genes" , flush = True )        
  Genes = {}
  for lName in lMembers: 
    _ , lCaseId , lMethod , _ = lName.name.split( "/" )  

    if lCaseId in ATRXmut: Index = 0
    elif lCaseId in ATRXwt: Index = 1
    elif lCaseId in ATRXother: Index = 2
    else: continue
    
    if lMethod == "RNA-Seq":
      for line in utf8reader( src.extractfile( lName ) ):
        line = [ i.strip() for i in line.split( "\t" ) ]
        if not line[0].startswith( "ENSG" ): continue
        lGeneName , lValue = line[1] , float( line[6] ) # Gene , TPM_unstranded
        if not lGeneName in Genes: Genes[ lGeneName ] = [ [] , [] , [] ]
        Genes[ lGeneName ][ Index ].append( lValue )
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( "Analysing and outputting" , flush = True )  
  dest.write( f"Gene\tCompare ATRXmut-ATRXwt\t\tATRXmut-ATRXother\t\tATRXwt-ATRXother\t\n" )
  dest.write( f"\tt-score\tp-value\tt-score\tp-value\tt-score\tp-value\n" )
  for k,v in sorted( Genes.items() ):
    t01, p01 = ttest_ind( v[0] , v[1] , equal_var = False )
    t02, p02 = ttest_ind( v[0] , v[2] , equal_var = False )
    t12, p12 = ttest_ind( v[1] , v[2] , equal_var = False )
    dest.write( f"{k}\t{t01}\t{p01}\t{t02}\t{p02}\t{t12}\t{p12}\n" )

    v0 , v1 , v2 = '\t'.join( [ str(i) for i in v[0] ] ) , '\t'.join( [ str(i) for i in v[1] ] ) , '\t'.join( [ str(i) for i in v[2] ] )    
    dest2.write( f"{k}\tATRX mutation\t{v0}\n" )
    dest2.write( f"\tATRX wild-type\t{v1}\n" )
    dest2.write( f"\tATRX others\t{v2}\n" )
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( "Complete" , flush = True )  
  # ----------------------------------------------------------------------------------------------------
