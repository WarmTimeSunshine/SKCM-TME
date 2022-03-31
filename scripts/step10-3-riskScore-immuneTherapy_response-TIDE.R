#Analyze the correlation of RiskScore with immunotherapy response
#（Predicted immunotherapy is from TIDE http://tide.dfci.harvard.edu./）

rm(list=ls())
options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(ggplot2)
library(ggpubr)
request_cancer = 'SKCM'

#step 01:To get the expression data to input into TIDE
if(T){
  workflow_type = "HTSeq - FPKM"
  exprSym_f=file.path(Rdata_dir,
                      paste(paste("TCGA",request_cancer,workflow_type,'exprSym_meta_groupT-N.RData',sep="-")))
  load(file = exprSym_f)
  exprSym_FPKM = na.omit(exprSym_FPKM)
  dim(exprSym_FPKM)
  dim(meta) 
  table(group_list_TN)
  
  #exprSym_FPKM procesing
  exprDat = exprSym_FPKM[,group_list_TN=="Tumor"]
  dim(exprDat)
  exprDat = as.data.frame(exprDat)
  exprDat_scale = as.data.frame(t(scale(t(exprDat))))
  write.table(exprDat_scale, file=file.path(Table_dir,'TIDE_SKCM_expr.txt'),
              sep = '\t', 
              quote = F,
              row.names=T, col.names=T)
}

#step 02:Boxplot showing the association of RS with responder and non-responder
if(T){
  phe_TCGA = read.csv(file = file.path(Table_dir,'TCGA_response_TIDE.csv'),
                            header = T)
  #Calculate RS for each patient
  l= load(file = file.path(Rdata_dir,'step50_LASSO_modeleGenes_coef.RData'))
  l
  genes = rownames(coef_out3)
  expr_lassoGenes_TCGA = log2(exprDat[genes,]+1)
  expr_lassoGenes_TCGA = as.data.frame(t(expr_lassoGenes_TCGA))
  expr_lassoGenes_TCGA$riskScore = coef_out3[1,1]*expr_lassoGenes_TCGA[,1] + 
    coef_out3[2,1]*expr_lassoGenes_TCGA[,2] +
    coef_out3[3,1]*expr_lassoGenes_TCGA[,3] + 
    coef_out3[4,1]*expr_lassoGenes_TCGA[,4] + 
    coef_out3[5,1]*expr_lassoGenes_TCGA[,5] + 
    coef_out3[6,1]*expr_lassoGenes_TCGA[,6] +
    coef_out3[7,1]*expr_lassoGenes_TCGA[,7]
  expr_lassoGenes_TCGA$Patient = rownames(expr_lassoGenes_TCGA)
  data_TCGA = merge(phe_TCGA[,c('Patient','Responder')],expr_lassoGenes_TCGA,by='Patient')
  data_TCGA$Response = ifelse(data_TCGA$Responder == 'True','R','NR')
    
  #t test and wilcox test
  x = data_TCGA$riskScore[data_TCGA$Responder == 'True']
  y = data_TCGA$riskScore[data_TCGA$Responder == 'False']
  wilcox.test(x,y,alternative="less",exact=TRUE,correct=TRUE)
  t.test(x,y,alternative="less")

  #Boxplot
  if(T){
    p <- ggboxplot(data_TCGA, 
                   x = 'Response', 
                   y = 'riskScore',
                   fill = 'Response', 
                   xlab = 'Predicted Response', 
                   ylab = 'Risk Score',
                   width = 0.5,
                   bxp.errorbar = T, 
                   bxp.errorbar.width = 0.2, 
                   outlier.shape = NA, 
                   order = c('R','NR'), 
                   palette = 'Dark2',
                   legend = "none" 
    )+ ylim(-2, 0.37)
    #Add significance markers
    my_comparisons <- list(c('NR','R'))
    p <- p + geom_signif(comparisons = my_comparisons,
                         test = "t.test",
                         vjust=0.4,#move the text up or down relative to the bracket
                         map_signif_level=c("*"=0.05,'ns'=1)
                         )
    p <- p + theme(
      plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
      plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
      plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
      axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
      axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
      axis.title.x  = element_text(color = 'black', size = 10, angle = 0),
      axis.title.y  = element_text(color = 'black', size = 10, angle = 90),
      legend.title  = element_blank(),
      legend.text   = element_text(color = 'black', size   = 7.5),
      axis.line = element_blank(), 
      panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
    )
    print(p)
    ggsave(filename = file.path(Figure_dir,paste0('boxplot_riskScore_response_TIDE-TCGA','_6-9cm.pdf')),
           plot = last_plot(), width = 6, height = 9, units = 'cm', dpi = 600)
  }
}













