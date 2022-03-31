#Survival analysis

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(survival)
library(survminer)

#load genes hubGenes_logFC_GS_MM
load(file = paste(Rdata_dir,"hubGenes_DEG-WGCNA.RData"))
#load clinical data (meta) and expression data (FPKM)
workflow_type="HTSeq - FPKM"
exprSym_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,workflow_type,'exprSym_meta_groupT-N.RData',sep="-")))
load(file = exprSym_f)
#3 files: exprSym_FPKM, meta, group_list_TN
dim(exprSym_FPKM) 

#Processing of dataframes meta and exprSym_FPKM
if(T){
  #exprSym_FPKM
  exprDat = exprSym_FPKM[rownames(hubGenes_logFC_GS_MM),group_list_TN=="Tumor"]
  exprDat = exprDat[,!duplicated(substr(colnames(exprDat),1,12))]
  dim(exprDat)
  exprDat = as.data.frame(exprDat)
  
  #meta：reserve ID,vital_status,days_to_death,days_to_last_follow_up of tumor patients
  sampleID = substr(colnames(exprDat),1,12)
  table(table(sampleID))
  clin_info = meta[,c('submitter_id','vital_status','days_to_death','days_to_last_follow_up')]
  clin_info = as.data.frame(clin_info)
  clin_info$time = ifelse(clin_info$vital_status=="Alive",
                          clin_info$days_to_last_follow_up, clin_info$days_to_death)
  clin_info$status=ifelse(clin_info$vital_status=="Dead",1,0)
  clin_info = clin_info[clin_info$submitter_id %in% sampleID,]
  phe = na.omit(clin_info[,c('submitter_id','vital_status','status','time')])
  phe = phe[phe$time>30,]

  exprSet = exprDat[,match(phe$submitter_id,substr(colnames(exprDat),1,12))]
  dim(exprSet)
}

t_pvalue = 0.05
#Survival analysis（Cox univariate regression）
if(T){
  colnames(phe)
  cox_results <-apply(exprSet , 1 , function(gene){
    phe$group=ifelse(gene>median(gene),2,1) #2=='high',1=='low'
    m=coxph(Surv(time, status) ~ group, data=phe)
    beta <- coef(m) 
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se
    tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
    return(tmp['group',])
  })
  cox_results=t(cox_results)
  cox_results = cox_results[order(cox_results[,4],decreasing = F),]
  write.csv(x=cox_results[,c(4,5,9,10)],
            file=file.path(Table_dir,
                           paste(paste("TCGA",request_cancer,'Uni_Cox_hubGenes.csv',sep="-"))),
            row.names = T, col.names=T)
  table(cox_results[,4]<t_pvalue) 
  sigGenes_cox = cox_results[cox_results[,4]<t_pvalue,]
  surv_f=file.path(Rdata_dir,
                   paste(paste("TCGA",request_cancer,'survival_results_cox.RData',sep="-")))
  save(cox_results,sigGenes_cox,file = surv_f)
  
  sigGenes_logFC_GS_MM = hubGenes_logFC_GS_MM[rownames(sigGenes_cox),]
}


