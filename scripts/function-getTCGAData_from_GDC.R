#Download mRNA and clinical information form TCGA
#For example (liver cancer):getTCGAData_from_GDC("LIHC","Transcriptome Profiling",
  # "Gene Expression Quantification", workflow_type="HTSeq - Counts",legacy=FALSE )

getTCGAData_from_GDC <- function(request_cancer,data_category,data_type, workflow_type,legacy){

  options(stringsAsFactors = F)
  Rdata_dir='../Rdata/'
  Figure_dir='../figures/'
  library(TCGAbiolinks)
  library(dplyr)
  #library(DT)
  library(SummarizedExperiment)
  
  #Download mRNA and clinical information
  tmpfile = file.path(Rdata_dir,
                      paste(paste("TCGA",request_cancer,data_type,workflow_type,'getTCGAData',sep="-"),
                            'RData',sep="."))
  if(!file.exists(tmpfile)){
    cancer_type <- paste("TCGA",request_cancer,sep="-")  
    print(cancer_type)
    #mRNA data
    query_TranscriptomeCounts <- GDCquery(project = cancer_type, 
                                          data.category = data_category, 
                                          data.type =  data_type, 
                                          workflow.type = workflow_type) 
    #"files.per.chunk" will make the API method only download n files at a time. 
    #This may reduce the download problems when the data size is too large
    GDCdownload(query_TranscriptomeCounts, method = "api", files.per.chunk = 10)  
    expdat <- GDCprepare(query = query_TranscriptomeCounts) 
    expr=assay(expdat)#expr is the expression data,rownames--ensemble ID,colnames--sample ID
    
    #clinical information
    meta <- GDCquery_clinic(project = paste("TCGA",request_cancer,sep="-") , type = "clinical")
    save(expr,meta,file = tmpfile)
  }
  
  #load mRNA and clinical data----
  load(file = tmpfile)
  print(dim(expr))
  print(dim(meta))
}

