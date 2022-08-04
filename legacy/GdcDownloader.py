from GdcLib import *

lFilter = Or(
            And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "RNA-Seq" ), Eq( "files.analysis.workflow_type" , "STAR - Counts"), Eq( "files.data_type" , "Gene Expression Quantification") ),
            And( Eq( "cases.project.project_id" , "TCGA-LGG" ), Eq( "files.experimental_strategy" , "WXS" ), Eq( "files.access" , "open") )
         )

lFileList = GetFileList( lFilter )
lTarball  = GetFiles( list( lFileList.keys() ) ) # Can add second argument (such as f"Raw-GDC-{TimeStamp}.tgz") to store raw
lFileName = SaveReformattedTarball( lFileList , lTarball )

           
