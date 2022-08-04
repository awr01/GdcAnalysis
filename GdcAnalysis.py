import argparse , tarfile , codecs , gzip
from scipy.stats import ttest_ind

utf8reader = codecs.getreader( 'utf-8' )

# Read the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')
args = parser.parse_args()


ATRXmut , ATRXwt , ATRXother = [] , [] , [] # Initialize empty lists to store case-ids for the different mutation-classifications
Genes = {}                                  # Initialize empty dictionary of gene-name to tpm-unstranded data per mutation-classification

print( f"Opening formatted tarball '{args.src}'" , flush=True )
with tarfile.open( args.src ) as src , open( args.dest , "w" ) as dest , open( "raw-" + args.dest , "w" ) as dest2: # Open the formatted tarball file and the destination files
  lMembers = src.getmembers()                               # Get the list of files in the tarball
  
  # ----------------------------------------------------------------------------------------------------  
  print( "Getting mutations" , flush = True )
  for lName in lMembers:                                    # Iterate over the files in the tarball
    _ , lCaseId , lMethod , _ = lName.name.split( "/" )     # Extract the case-id and file-type from the name     
    if lMethod == "WXS":                                    # Only look at WXS files
      print( f"." , end="" , flush=True )        
      Mutations = []                                        # For each file initialize an empty mutations list
      
      for line in utf8reader( gzip.GzipFile( fileobj = src.extractfile( lName ) ) ): # Iterate over each line in the file
        if line[0] == "#" or line.startswith( "Hugo_Symbol" ) : continue             # Ignore comments and headers
        line = [ i.strip() for i in line.split( "\t" , maxsplit = 9 ) ]              # Split the line at tabs up to where we need it
        if line[0] == "ATRX" : Mutations.append( line[8] )                           # If the mutation is on ATRX, add it to the list

      if len( Mutations ) == 0:                             # If there are no mutations
        ATRXwt.append( lCaseId )                            # add case-id to the wild-type list
        continue                                            # and go to the next file
      
      Mutations = list( filter( lambda a: (a != "Silent") and (a!="Splice_Site") , Mutations ) ) # Remove all silent and splice-site from the mutation list
      if len( Mutations ) == 0:                             # If there are now no mutations (i.e. there were only silents or slices)
        ATRXother.append( lCaseId )                         # add case-id to the "others" list
        continue                                            # and go to the next file

      ATRXmut.append( lCaseId )                             # Else add case-id to the mutations list
  print( flush = True ) 
  # ----------------------------------------------------------------------------------------------------
     
  # ----------------------------------------------------------------------------------------------------
  print( "Getting genes" , flush = True )        
  for lName in lMembers:                                    # Iterate over the files in the tarball
    _ , lCaseId , lMethod , _ = lName.name.split( "/" )     # Extract the case-id and file-type from the name  

    if lMethod == "RNA-Seq":                                # Only look at RNA-Seq files
      print( f"." , end="" , flush=True )        

      if   lCaseId in ATRXmut:   Index = 0                  # Classify the case-id by its mutation-classification      
      elif lCaseId in ATRXwt:    Index = 1
      elif lCaseId in ATRXother: Index = 2
      else: continue                                        # If there is no WXS data for this case-id, move to the next case
    
      for line in utf8reader( src.extractfile( lName ) ):   # Iterate over each line in the file
        line = [ i.strip() for i in line.split( "\t" , maxsplit = 7 ) ] # Split the line at tabs up to where we need it
        if not line[0].startswith( "ENSG" ): continue       # Ignore comments and headers
        lGeneName , lValue = line[1] , float( line[6] )     # Store the gene-name and TPM_unstranded data as variables
        if not lGeneName in Genes: Genes[ lGeneName ] = [ [] , [] , [] ] # If we haven't seen this gene before, initialize new entry in the dictionary
        Genes[ lGeneName ][ Index ].append( lValue )        # Append the tpm-unstranded data to the list corresponding to this gene and this ATRX mutation-classification
  print( flush = True ) 
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( "Analysing and outputting" , flush = True )  
  dest.write( f"Gene\tCompare ATRXmut-ATRXwt\t\tATRXmut-ATRXother\t\tATRXwt-ATRXother\t\n" ) # Write headers
  dest.write( f"\tt-score\tp-value\tt-score\tp-value\tt-score\tp-value\n" )
  
  for k,v in sorted( Genes.items() ):                       # For each gene
    t01, p01 = ttest_ind( v[0] , v[1] , equal_var = False ) # Calculate the t-score between each pair of lists
    t02, p02 = ttest_ind( v[0] , v[2] , equal_var = False )
    t12, p12 = ttest_ind( v[1] , v[2] , equal_var = False )
    dest.write( f"{k}\t{t01}\t{p01}\t{t02}\t{p02}\t{t12}\t{p12}\n" ) # Write t-score and p-value data to file

    # Copy the raw-data to tab-delimeted strings
    v0 , v1 , v2 = '\t'.join( [ str(i) for i in v[0] ] ) , '\t'.join( [ str(i) for i in v[1] ] ) , '\t'.join( [ str(i) for i in v[2] ] )    
    dest2.write( f"{k}\tATRX mutation\t{v0}\n" )            # Write the raw data to "raw" file
    dest2.write( f"\tATRX wild-type\t{v1}\n" )
    dest2.write( f"\tATRX others\t{v2}\n" )
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( "Complete" , flush = True )  
  # ----------------------------------------------------------------------------------------------------
