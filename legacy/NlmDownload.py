import requests, json, os , bz2 , _pickle


def GeneDetails( aGeneList ):
  if isinstance( aGeneList , str ) : aGeneList = [ aGeneList ]

  Filename = "genes.nlm.nih.pkl.bz2"
  Data = {}

  if os.path.isfile( Filename ):
    with bz2.BZ2File( Filename , 'r' ) as src: Data = _pickle.load( src )  
  
  lRequired = [ f'"{i}"' for i in aGeneList if not i in Data ]
 
  if len(lRequired):
    print( f"Downloading info for {len(lRequired)} genes" , flush = True )
    for i in range( 0 , len( lRequired ) , 1000 ):
      print( f"." , end="" , flush=True )        
      lRequest = ','.join( lRequired[ i : min( len( lRequired ) , i+1000 ) ] )
      while True:
        lResponse = requests.post( "https://api.ncbi.nlm.nih.gov/datasets/v1/gene" , data = '{"symbols_for_taxon":{"symbols":[ '+lRequest+' ],"taxon":"Homo sapiens"}}' , headers = { "Accept" : "application/json" , "Content-Type" : "application/json" } )     
        if lResponse.status_code == 200 : break
        print( f"!" , end="" , flush=True )        
        
      lResponse = json.loads( lResponse.content.decode("utf-8") )
      if "genes" in lResponse:
        for i in lResponse["genes"]: 
          if "gene" in i:
            i["gene"].pop('transcripts', None)
            Data[ i["query"][0] ] = i["gene"]
          else: 
            Data[ i["query"][0] ] = None

    for i in aGeneList:
      if not i in Data: Data[ i ] = None

    print( flush=True )
    with bz2.BZ2File( Filename , 'w' ) as dest: _pickle.dump( Data , dest )

  return { i:Data[i] for i in aGeneList }


# GeneDetails( [ "ATRX" , "GAPDH" , "TUBA1A" , "ACTB" , "ALB" , "ALOX12" , "ANGPTL7" , "AOX1" , "APOE" , "ATOX1" , 
                    # "BNIP3" , "CAT" , "CCL5" , "CCS" , "CSDE1" , "CYBA" , "CYGB" , "DGKK" , "DHCR24" , "DUOX1" , "DUOX2" , 
                    # "DUSP1" , "EPHX2" , "EPX" , "FOXM1" , "GLRX2" , "GPR156" , "GPX1" , "GPX2" , "GPX3" , "GPX4" , "GPX5" , 
                    # "GPX6" , "GPX7" , "GSR" , "GSS" , "GSTZ1" , "GTF2I" , "KRT1" , "LPO" , "MBL2" , "MGST3" , "MPO" , 
                    # "MPV17" , "MSRA" , "MT3" , "TESMIN" , "NCF1" , "NCF1B" , "NCF1C" , "NCF2" , "NME5" , "NOS2" , "NOX5" , 
                    # "NUDT1" , "OXR1" , "OXSR1" , "PDLIM1" , "IPCEF1" , "PNKP" , "PRDX1" , "PRDX2" , "PRDX3" , "PRDX4" , "PRDX5" , "PRDX6" , 
                    # "PREX1" , "PRG3" , "PRNP" , "PTGS1" , "PTGS2" , "PXDN" , "PXDNL" , "RNF7" , "SCARA3" , "SELENOS" , "SELENOP" , "SFTPD" , 
                    # "SGK2" , "SIRT2" , "SOD1" , "SOD2" , "SOD3" , "SRXN1" , "STK25" , "TPO" , "TTN" , "TXNDC2" , "TXNRD1" , "TXNRD2" ] )
                    
# for i,j in GetGenes( "ATRX" ).items() :
  # for k,l in j.items() :
    # print( k , l )