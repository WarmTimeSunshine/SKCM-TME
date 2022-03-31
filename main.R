
# Set work disc
if(T){
  setwd('scripts/')
}

#step00: Install packages
if(F){
  source('step00-install-packages.R')
}

# STEP 01:TCGA download
if(F){
  rm(list=ls())
  data_category="Transcriptome Profiling"
  data_type="Gene Expression Quantification" 
  workflow_type="HTSeq - Counts"
  source('function-getTCGAData_from_GDC.R')
  request_cancer="SKCM"
  getTCGAData_from_GDC(request_cancer,data_category,
                       data_type, workflow_type,legacy=FALSE)

  #gene ID convert（ensemble ID---symbol）
  source("step01-2-IDTrans-by-bioMart.R")
  source("step01-3-checkData.R")
}

##STEP 02:ESTIMATE Score
if(T){
  rm(list = ls())  
  request_cancer="SKCM"
  source('step02-1-ESTIMATEScore.R')

  rm(list=ls())
  request_cancer="SKCM"
  source('step02-2-boxplot_ESTIMATEScore-stage.R')
  rm(list=ls())
  request_cancer="SKCM"
  source('step02-3-surv_ESTIMATEScore.R')
  rm(list=ls())
  request_cancer="SKCM"
  source('step02-4-boxplot_ESTIMATEScore-mut.R')
}

#STEP 03:DEG:Differentially Expressed Genes
if(T){
  rm(list = ls())  
  request_cancer="SKCM"
  source('step03-DEG-score_high-low.R')
}

#STEP 05:WGCNA，to identify immune-related modules
if(T){
  # re-download expression data (FPKM) for WGCNA
  if(T){
    rm(list=ls())
    request_cancer="SKCM"
    data_category="Transcriptome Profiling"
    data_type="Gene Expression Quantification"
    workflow_type="HTSeq - FPKM"
    source('function-getTCGAData_from_GDC.R')
    getTCGAData_from_GDC(request_cancer,data_category,
                         data_type, workflow_type,legacy=FALSE)

    source("step05-0-IDTrans-by-bioMart.R")
    #expression data called exprSym_FPKM,and clinical information called meta,
    #group information called group_list_TN（Tumor,Normal）
    #all saved as'TCGA-SKCM-FPKM_exprSym_meta_groupT-N.Rdata'
  }
  
  #STEP 05: WGCNA
  rm(list=ls())
  request_cancer="SKCM"
  source('step05-WGCNA.R')
  
  #STEP 06:combination of DEGs and WGCNA to identify immune-related genes
  rm(list=ls())
  request_cancer="SKCM"
  source('step06-commonGenes_DEG-WGCNA.R')
  
  #STEP 08:enrichment analysis（GO、KEGG）
  rm(list=ls())
  request_cancer="SKCM"
  source('step08-anno-go-kegg-GSEA.R')
}

#STEP 07:prognostic value of screened genes
if(T){
  rm(list=ls())
  request_cancer="SKCM"
  source('step07-survivalAnalysis.R')
}

#STEP 09:hub genes selection；LASSO model，ROC curve
if(T){
  #hub genes selection
  rm(list=ls())
  request_cancer="SKCM"
  source('step09-0-hubGenes_from_TCGA-GEO.R')
  
  #LASSO
  rm(list=ls())
  request_cancer="SKCM"
  source('step09-1-LASSOmodele-construction_hubGenes_from_TCGA-GEO.R')
  
  # TCGA：KM curve
  rm(list=ls())
  request_cancer="SKCM"
  source('step09-2-LASSOmodele-KM.R')
  
  # TCGA：ROC curve
  rm(list=ls())
  request_cancer="SKCM"
  source('step09-3-LASSOmodele-ROC.R')
  
  # GEO：KM curve
  rm(list=ls())
  GSEName='GSE65904'
  source('step09-4-LASSOmodele-KM_GEO.R')
  
  #GEO：ROC curve
  rm(list=ls())
  GSEName = "GSE65904"
  source('step09-5-LASSOmodele-ROC_GEO.R')
  
  #Cox analysis; nomogram
  rm(list=ls())
  request_cancer="SKCM"
  source('step09-6-LASSOmodele-MultiCox_nomogram.R')
}

#step 10:Association of RS with immunotherapy
if(T){
  rm(list=ls())
  request_cancer="SKCM"
  source('step10-1-LASSOmodele-immuneInfiltration.R')
  
  rm(list=ls())
  request_cancer="SKCM"
  source('step10-2-LASSOmodele-ICIgenes_expr.R')
  
  #predict immunothepary response from TIDE website of TCGA-SKCM
  source('step10-3-riskScore-immuneTherapy_response-TIDE.R')
}



