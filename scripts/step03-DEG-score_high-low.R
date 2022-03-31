#DEGs between Immune Score-high and -low groups(DESeq2)

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'

# To get estimate_score
estimate_score = read.table(file = file.path(Table_dir, paste("ESTIMATE",request_cancer,'fromWeb.txt',sep="_")),
                            header=TRUE)
dim(estimate_score) 

exprSym_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,'exprSym_meta_groupT-N.RData',sep="-")))
load(file = exprSym_f)
table(group_list_TN)
dim(exprSym)

exprSet = na.omit(exprSym[,group_list_TN=="Tumor"])
dim(exprSet)
colnames(exprSet) = substr(colnames(exprSet),1,15)
table(table(colnames(exprSet))) 

samples = intersect(colnames(exprSet),estimate_score$ID)
exprSet = exprSet[,samples]
dim(exprSet)
estimate_score = estimate_score[estimate_score$ID %in% samples,]
dim(estimate_score)

# Divide patients into 2 groups based on median of immune score
if(T){
  groupBy = "Immune_score"
  group_list_ImmuneScore=ifelse(estimate_score[,groupBy]>
                                  median(estimate_score[,groupBy]),
                                'High','Low')
  names(group_list_ImmuneScore) = estimate_score$ID
  table(group_list_ImmuneScore)
}

#DEG based on group_list_ImmuneScore----
if(T){
  library(DESeq2)
  tmp_f=file.path(Rdata_dir,'DEG_by_ImmuneScore.Rdata')
  if(!file.exists(tmp_f)){
    identical(colnames(exprSet),names(group_list_ImmuneScore))
    if(!identical(colnames(exprSet),names(group_list_ImmuneScore))){
      exprSet = exprSet[,names(group_list_ImmuneScore)]
    }
    
    #To get colData
    colData <- data.frame(row.names=colnames(exprSet), 
                          group_list=group_list_ImmuneScore) 
    exprSet_interger <- round(exprSet)
    rownames(exprSet_interger) = rownames(exprSet)
    colnames(exprSet_interger) = colnames(exprSet)
    dds <- DESeqDataSetFromMatrix(countData = exprSet_interger,
                                  colData = colData,
                                  design = ~ group_list)
    dds <- DESeq(dds)

    res <- results(dds, contrast=c("group_list","High","Low"))
    resOrdered <- res[order(res$padj),]
    head(resOrdered)
    DEG =as.data.frame(resOrdered)
    DESeq2_DEG = na.omit(DEG)
    DESeq2_DEG_by_ImmuneScore=DESeq2_DEG[,c(2,6)]
    colnames(DESeq2_DEG_by_ImmuneScore)=c('log2FoldChange','pvalue')  
    save(DESeq2_DEG_by_ImmuneScore,file = tmp_f)
  }
  load(tmp_f)
}

logFC_t=1
pvalue_t=0.05
if(T){
  need_DEG_by_ImmuneScore = DESeq2_DEG_by_ImmuneScore[
    DESeq2_DEG_by_ImmuneScore$pvalue<pvalue_t & 
      (abs(DESeq2_DEG_by_ImmuneScore$log2FoldChange)>logFC_t),]
  UP_DEG_by_ImmuneScore = need_DEG_by_ImmuneScore[need_DEG_by_ImmuneScore$log2FoldChange>logFC_t,]
  DOWN_DEG_by_ImmuneScore = need_DEG_by_ImmuneScore[need_DEG_by_ImmuneScore$log2FoldChange<(0-logFC_t),]
  dim(need_DEG_by_ImmuneScore)
  dim(UP_DEG_by_ImmuneScore)
  dim(DOWN_DEG_by_ImmuneScore)
  save(need_DEG_by_ImmuneScore, UP_DEG_by_ImmuneScore,DOWN_DEG_by_ImmuneScore, 
       file = file.path(Rdata_dir,'need_DEG_by_ImmuneScore.Rdata'))
  
  #volcano
  load(file = file.path(Rdata_dir,'DEG_by_ImmuneScore.Rdata'))
  source('function-draw_volcano.R')
  draw_volcano(DEG_result=DESeq2_DEG_by_ImmuneScore,
               pvalue_t=pvalue_t,logFC_t=logFC_t,method='DESeq2',figure_dir=Figure_dir)
}



