#Survival analysis of hubgenes based on TCGA and GEO datasets

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(survival)
library(survminer)

#load genes hubGenes_logFC_GS_MM
load(file = paste(Rdata_dir,"hubGenes_DEG-WGCNA.RData"))
workflow_type="HTSeq - FPKM"
exprSym_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,workflow_type,'exprSym_meta_groupT-N.RData',sep="-")))
load(file = exprSym_f)
#3 files included：exprSym_FPKM, meta, group_list_TN
dim(exprSym_FPKM) 

#processing of datdaframes meta and exprSym_FPKM
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
  table(cox_results[,4]<t_pvalue)
  sigGenes_cox_TCGA = cox_results[cox_results[,4]<t_pvalue,]
  surv_f=file.path(Rdata_dir,
                   paste(paste("TCGA",request_cancer,'survival_results_cox.RData',sep="-")))
  save(cox_results,sigGenes_cox_TCGA,file = surv_f)
  
  sigGenes_logFC_GS_MM = hubGenes_logFC_GS_MM[rownames(sigGenes_cox_TCGA),]
}

#Survival analysis based on GEO datadset
if(T){
  #Procesing of meta and exprSet
  if(T){
    GSEName = 'GSE65904'
    RData_file=file.path(Rdata_dir,paste0(GSEName,'.RData'))
    load(file = RData_file)
    exprDat = exprSet[rownames(exprSet) %in% rownames(sigGenes_cox_TCGA),]
    dim(exprDat)
    exprDat = as.data.frame(exprDat)
    
    clin_info = meta[,c(2,3)]
    clin_info = as.data.frame(clin_info)
    colnames(clin_info) = c('status','time')
    phe = clin_info[clin_info$time>30,]

    sample_ID = intersect(rownames(phe),colnames(exprDat))
    exprSet = exprDat[,sample_ID]
    phe = phe[sample_ID,]
    dim(exprSet)
    dim(phe)
    identical(rownames(phe),colnames(exprSet))
  }
  
  #Survival analysis（Cox univariate regression）
  if(T){
    colnames(phe)
    cox_results <-apply(exprSet , 1 , function(gene){
      phe$group=ifelse(gene>median(gene),2,1) 
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
    table(cox_results[,4]<t_pvalue)
    sigGenes_cox_GEO = cox_results[cox_results[,4]<t_pvalue,]
    sigGenes_cox_TCGA = sigGenes_cox_TCGA[rownames(sigGenes_cox_GEO),]
    write.csv(sigGenes_cox_TCGA[,c(4,5,9,10)],
              file=file.path(Table_dir,"sigGenes_HR_p_TCGA.csv"),
              row.names = T, col.names=T)
    write.csv(sigGenes_cox_GEO[,c(4,5,9,10)],
              file=file.path(Table_dir,paste0("sigGenes_HR_p_GEO_",GSEName,".csv")),
              row.names = T, col.names=T)
    surv_f=file.path(Rdata_dir,
                     paste('survival_results_cox_TCGA.RData',GSEName,sep="-"))
    save(sigGenes_cox_TCGA,sigGenes_cox_GEO,file = surv_f)
    
    sigGenes_logFC_GS_MM = hubGenes_logFC_GS_MM[rownames(sigGenes_cox_GEO),]
  }
}










