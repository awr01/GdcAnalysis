import argparse , tarfile , codecs , gzip
from scipy.stats import ttest_ind
from numpy import mean , var

utf8reader = codecs.getreader( 'utf-8' )

# Read the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--src' , required=True , help='The source tarball')
parser.add_argument( '--dest' , required=True , help='The destination file')
args = parser.parse_args()


GenesOfInterest = [ "ATRX" , "GAPDH" , "TUBA1A" , "ACTB" , "ALB" , "ALOX12" , "ANGPTL7" , "AOX1" , "APOE" , "ATOX1" , 
                    "BNIP3" , "CAT" , "CCL5" , "CCS" , "CSDE1" , "CYBA" , "CYGB" , "DGKK" , "DHCR24" , "DUOX1" , "DUOX2" , 
                    "DUSP1" , "EPHX2" , "EPX" , "FOXM1" , "GLRX2" , "GPR156" , "GPX1" , "GPX2" , "GPX3" , "GPX4" , "GPX5" , 
                    "GPX6" , "GPX7" , "GSR" , "GSS" , "GSTZ1" , "GTF2I" , "KRT1" , "LPO" , "MBL2" , "MGST3" , "MPO" , 
                    "MPV17" , "MSRA" , "MT3" , "TESMIN" , "NCF1" , "NCF1B" , "NCF1C" , "NCF2" , "NME5" , "NOS2" , "NOX5" , 
                    "NUDT1" , "OXR1" , "OXSR1" , "PDLIM1" , "IPCEF1" , "PNKP" , "PRDX1" , "PRDX2" , "PRDX3" , "PRDX4" , "PRDX5" , "PRDX6" , 
                    "PREX1" , "PRG3" , "PRNP" , "PTGS1" , "PTGS2" , "PXDN" , "PXDNL" , "RNF7" , "SCARA3" , "SELENOS" , "SELENOP" , "SFTPD" , 
                    "SGK2" , "SIRT2" , "SOD1" , "SOD2" , "SOD3" , "SRXN1" , "STK25" , "TPO" , "TTN" , "TXNDC2" , "TXNRD1" , "TXNRD2" ]



ATRXmut , ATRXwt = [] , []  # Initialize empty lists to store case-ids for the different mutation-classifications
Genes = {}                  # Initialize empty dictionary of gene-name to tpm-unstranded data per mutation-classification

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
        continue                                            # go to the next file

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
      else: continue                                        # If there is no WXS data for this case-id, move to the next case
    
      for line in utf8reader( src.extractfile( lName ) ):   # Iterate over each line in the file
        line = [ i.strip() for i in line.split( "\t" , maxsplit = 7 ) ] # Split the line at tabs up to where we need it
        if not line[0].startswith( "ENSG" ): continue       # Ignore comments and headers
        lGeneName , lValue = line[1] , float( line[6] )     # Store the gene-name and TPM_unstranded data as variables
        if not lGeneName in Genes: Genes[ lGeneName ] = [ [] , [] ] # If we haven't seen this gene before, initialize new entry in the dictionary
        Genes[ lGeneName ][ Index ].append( lValue )        # Append the tpm-unstranded data to the list corresponding to this gene and this ATRX mutation-classification
  print( flush = True ) 
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( "Analysing and outputting" , flush = True )  
  dest.write( f"Gene\tFlagged\tMut count\tMut mean\tMut var\tWT count\tWT mean\tWT var\t Mut mean/WT mean\tATRXmut:ATRXwt t-score\tATRXmut:ATRXwt p-value\n" ) # Write headers
  
  for k,v in sorted( Genes.items() ):                       # For each gene 
    num0 , mean0 , var0 = len( v[0] ) , mean( v[0] ) , var( v[0] )               # Calculate the mean and variance of the muts
    num1 , mean1 , var1 = len( v[1] ) , mean( v[1] ) , var( v[1] )               # Calculate the mean and variance of the WTs
    
    if mean0 == 0 or mean1 == 0 : continue
    
    t01, p01 = ttest_ind( v[0] , v[1] , equal_var = False ) # Calculate the t-score between each pair of lists
    meanratio = mean0 / mean1

    flag = "*" if k in GenesOfInterest else ""
      
    
    dest.write( f"{k}\t{flag}\t{num0}\t{mean0}\t{var0}\t{num1}\t{mean1}\t{var1}\t{meanratio}\t{t01}\t{p01}\n" ) # Write everything to file

    # Copy the raw-data to tab-delimeted strings
    v0 , v1 = '\t'.join( [ str(i) for i in v[0] ] ) , '\t'.join( [ str(i) for i in v[1] ] )
    dest2.write( f"{k}\t{flag}\tATRX mutation\t{v0}\n" )            # Write the raw data to "raw" file
    dest2.write( f"\t\tATRX wild-type\t{v1}\n" )
  # ----------------------------------------------------------------------------------------------------

  # ----------------------------------------------------------------------------------------------------
  print( "Complete" , flush = True )  
  # ----------------------------------------------------------------------------------------------------
