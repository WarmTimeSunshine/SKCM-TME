#Analyze correlation between riskScore and ICI genes(PDCD1/PD1, CD274/PD-L1 and CTLA-4)
#expression in TCGA and GEO datasets

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(ggplot2)
library(ggpubr)

#1.TCGA dataset
#1.1 data procesing
if(T){
  workflow_type = "HTSeq - FPKM"
  exprSym_f=file.path(Rdata_dir,
                      paste(paste("TCGA",request_cancer,workflow_type,'exprSym_meta_groupT-N.RData',sep="-")))
  load(file = exprSym_f)
  exprSym_FPKM = na.omit(exprSym_FPKM)
  dim(exprSym_FPKM)
  dim(meta) 
  table(group_list_TN)
  
  #exprSym_FPKM
  exprDat = exprSym_FPKM[,group_list_TN=="Tumor"]
  dim(exprDat)
  exprDat = as.data.frame(exprDat)
  
  #Calculate riskScore for each patient
  load(file = file.path(Rdata_dir,'step50_LASSO_modeleGenes_coef.RData'))
  genes = rownames(coef_out3)
  expr_lassoGenes = log2(exprDat[genes,]+1)
  expr_lassoGenes = as.data.frame(t(expr_lassoGenes))
  expr_lassoGenes$riskScore = coef_out3[1,1]*expr_lassoGenes[,1] + 
    coef_out3[2,1]*expr_lassoGenes[,2] +
    coef_out3[3,1]*expr_lassoGenes[,3] + 
    coef_out3[4,1]*expr_lassoGenes[,4] + 
    coef_out3[5,1]*expr_lassoGenes[,5] + 
    coef_out3[6,1]*expr_lassoGenes[,6] +
    coef_out3[7,1]*expr_lassoGenes[,7]
  expr_lassoGenes$RS_group = ifelse(expr_lassoGenes$riskScore>median(expr_lassoGenes$riskScore),
                                    'High_risk','Low_risk')
  
  x = exprSym_FPKM[,group_list_TN=="Tumor"]
  expr_log = log2(x+1)
  x=t(expr_log[c('PDCD1','CD274','CTLA4'),])
  identical(rownames(x),rownames(expr_lassoGenes))
  expr = cbind(expr_lassoGenes, x)
  cor(expr$riskScore,expr[,c('PDCD1','CD274','CTLA4')],method = "pearson")
  expr_plot = expr[,8:12]
  
  #Collect data required to draw correlation graph
  if(T){
    sample_rep = rep(rownames(expr_plot),3)
    RS_rep = rep(expr_plot$riskScore,3)
    z = array(colnames(expr_plot[,-c(1,2)]),dim=c(ncol(expr_plot[,-c(1,2)]),1))
    gene_matrix = apply(z, 1, function(x) rep(x,nrow(expr_plot)))
    gene_rep = as.vector(gene_matrix)
    expr_dat = as.vector(as.matrix(expr_plot[,-c(1,2)]))
    
    cor_plot = as.data.frame(cbind(sample_rep,RS_rep,gene_rep,expr_dat))
    cor_plot$expr_dat = as.numeric(cor_plot$expr_dat)
    cor_plot$RS_rep = as.numeric(cor_plot$RS_rep)
    cor_plot$gene_rep = factor(cor_plot$gene_rep,levels = c('PDCD1','CD274','CTLA4'))
  }
  #Collect data required to draw boxplot
  if(T){
    sample_rep = rep(rownames(expr_plot),3)
    group_rep = rep(expr_plot$RS_group,3)
    z = array(colnames(expr_plot[,-c(1,2)]),dim=c(ncol(expr_plot[,-c(1,2)]),1))
    gene_matrix = apply(z, 1, function(x) rep(x,nrow(expr_plot)))
    gene_rep = as.vector(gene_matrix)
    expr_dat = as.vector(as.matrix(expr_plot[,-c(1,2)]))
    
    box_plot = as.data.frame(cbind(sample_rep,gene_rep,group_rep,expr_dat))
    box_plot$expr_dat = as.numeric(box_plot$expr_dat)
    box_plot$gene_rep = factor(box_plot$gene_rep,levels = c('PDCD1','CD274','CTLA4'))
    box_plot$group_rep = factor(box_plot$group_rep,levels = c('High_risk','Low_risk'))
  }
}
#1.2 correlation graph
if(T){
  library(ggpubr)
  p1 <- ggscatter(cor_plot,
                  x = "RS_rep", 
                  y = "expr_dat",
                  color = 'black', #color of points
                  size =0.2,  #size of points
                  add = "reg.line",
                  add.params = list(color = "blue",size=0.6),
                  conf.int = F,
                  cor.coef = TRUE, 
                  cor.method = "pearson",
                  cor.coef.size = 2,
                  xlab = "Risk Score", 
                  ylab = "Expression",
                  ggtheme =  theme(
                    panel.background = element_rect(fill = 'white'), 
                    legend.key = element_rect(fill = 'white'),
                    plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
                    plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
                    plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
                    axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
                    axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
                    axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
                    axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
                    legend.title  = element_text(color = 'black', size  = 6),
                    legend.text   = element_text(color = 'black', size   = 6),
                    axis.line = element_blank(), 
                    panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
                  )
  )+ylim(0.0,8)
  p <- p1+facet_wrap(.~gene_rep, 
                     scales = 'fixed',
                     nrow =1 
  )+theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text  = element_text(color = 'black',size = 7.5, angle = 0)
  )
  print(p)
  ggsave(filename = file.path(Figure_dir,'cor_riskScore-ICIgenes_expr-TCGA_9-3.5cm.pdf'),
         plot = last_plot(), width = 9, height = 3.5, 
         units = 'cm', dpi = 600)
}
#1.3 Boxplot
if(T){
  p1 <- ggboxplot(box_plot, 
                  x = 'gene_rep', 
                  y = 'expr_dat', 
                  fill = 'group_rep', 
                  xlab = FALSE,
                  ylab = 'Expression',
                  width = 0.6,#numeric value between 0 and 1 specifying box width.
                  bxp.errorbar = T, 
                  bxp.errorbar.width = 0.2, 
                  outlier.shape = NA, 
                  palette = 'npg' 
  ) + ylim(0, 6.5)
  p1 <- p1 + scale_x_discrete(
    breaks = c('PDCD1','CD274','CTLA4'),
    label = c('PD1','PD-L1','CTLA4'))
  p <- p1 + theme(
    plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 7.5, angle = 45,
                                 hjust = 0.5, vjust = 0.5),
    axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
    axis.title.x  = element_blank(),
    axis.title.y  = element_text(color = 'black', size = 10, angle = 90),
    legend.title  = element_blank(),
    legend.text   = element_text(color = 'black', size   = 7.5),
    legend.position ="top",
    axis.line = element_blank(), 
    panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
  )
  print(p)
  ggsave(filename = file.path(Figure_dir,'boxplot_riskScore_ICIgenes_expr_TCGA_4.5-9.3cm.pdf'),
         plot = p, width = 4.5, height = 9.3, units = 'cm', dpi = 600)
}

#2.GEO dataset---GSE65904
#2.1 data collection
if(T){
  GSEName = 'GSE65904'
  RData_file=file.path(Rdata_dir,paste0(GSEName,'.RData'))
  l = load(file = RData_file)
  l
  exprDat = as.data.frame(exprSet)
  dim(exprDat)
  
  #Calculate riskScore for each sample
  load(file = file.path(Rdata_dir,'step50_LASSO_modeleGenes_coef.RData'))
  genes = rownames(coef_out3)
  expr_lassoGenes = exprDat[genes,]
  expr_lassoGenes = as.data.frame(t(expr_lassoGenes))
  expr_lassoGenes$riskScore = coef_out3[1,1]*expr_lassoGenes[,1] + 
    coef_out3[2,1]*expr_lassoGenes[,2] +
    coef_out3[3,1]*expr_lassoGenes[,3] + 
    coef_out3[4,1]*expr_lassoGenes[,4] + 
    coef_out3[5,1]*expr_lassoGenes[,5] + 
    coef_out3[6,1]*expr_lassoGenes[,6] +
    coef_out3[7,1]*expr_lassoGenes[,7]
  expr_lassoGenes$RS_group = ifelse(expr_lassoGenes$riskScore>median(expr_lassoGenes$riskScore),
                                    'High_risk','Low_risk')
  
  x=t(exprDat[c('PDCD1','CD274','CTLA4'),])
  identical(rownames(x),rownames(expr_lassoGenes))
  expr = cbind(expr_lassoGenes, x)
  cor(expr$riskScore,expr[,c('PDCD1','CD274','CTLA4')],method = "pearson")
  expr_plot = expr[,8:12]
  
  if(T){
    sample_rep = rep(rownames(expr_plot),3)
    RS_rep = rep(expr_plot$riskScore,3)
    z = array(colnames(expr_plot[,-c(1,2)]),dim=c(ncol(expr_plot[,-c(1,2)]),1))#列名转为向量
    gene_matrix = apply(z, 1, function(x) rep(x,nrow(expr_plot)))
    gene_rep = as.vector(gene_matrix)
    expr_dat = as.vector(as.matrix(expr_plot[,-c(1,2)]))
    
    cor_plot = as.data.frame(cbind(sample_rep,RS_rep,gene_rep,expr_dat))
    cor_plot$expr_dat = as.numeric(cor_plot$expr_dat)
    cor_plot$RS_rep = as.numeric(cor_plot$RS_rep)
    cor_plot$gene_rep = ifelse(cor_plot$gene_rep=='PDCD1','PD1',
                               ifelse(cor_plot$gene_rep=='CD274','PD-L1','CTLA4'))
    cor_plot$gene_rep = factor(cor_plot$gene_rep,levels = c('PD1','PD-L1','CTLA4'))
  
    }
  if(T){
    sample_rep = rep(rownames(expr_plot),3)
    group_rep = rep(expr_plot$RS_group,3)
    z = array(colnames(expr_plot[,-c(1,2)]),dim=c(ncol(expr_plot[,-c(1,2)]),1))#列名转为向量
    gene_matrix = apply(z, 1, function(x) rep(x,nrow(expr_plot)))
    gene_rep = as.vector(gene_matrix)
    expr_dat = as.vector(as.matrix(expr_plot[,-c(1,2)]))
    
    box_plot = as.data.frame(cbind(sample_rep,gene_rep,group_rep,expr_dat))
    box_plot$expr_dat = as.numeric(box_plot$expr_dat)
    box_plot$gene_rep = factor(box_plot$gene_rep,levels = c('PDCD1','CD274','CTLA4'))
    box_plot$group_rep = factor(box_plot$group_rep,levels = c('High_risk','Low_risk'))
  }
}
#1.2 correlation plot
if(T){
  library(ggpubr)
  p1 <- ggscatter(cor_plot,
                  x = "RS_rep",
                  y = "expr_dat",
                  color = 'black',
                  size =0.2, 
                  add = "reg.line",
                  add.params = list(color = "blue",size=0.6),
                  conf.int = F,
                  cor.coef = TRUE, 
                  cor.method = "pearson",
                  cor.coef.size = 2,
                  xlab = "Risk Score",
                  ylab = "Expression",
                  ggtheme =  theme(
                    panel.background = element_rect(fill = 'white'), 
                    legend.key = element_rect(fill = 'white'),
                    plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
                    plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
                    plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
                    axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
                    axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
                    axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
                    axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
                    legend.title  = element_text(color = 'black', size  = 6),
                    legend.text   = element_text(color = 'black', size   = 6),
                    axis.line = element_blank(), 
                    panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
                  )
  )
  p <- p1+facet_wrap(.~gene_rep, 
                     scales = 'fixed',
                     nrow =1 
  )+theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text  = element_text(color = 'black',size = 7.5, angle = 0)
  )
  print(p)
  ggsave(filename = file.path(Figure_dir,'cor_riskScore-ICIgenes_expr-GSE65904_9-3.5cm.pdf'),
         plot = last_plot(), width = 9, height = 3.5, 
         units = 'cm', dpi = 600)
}
#1.3 boxplot 
if(T){
  p1 <- ggboxplot(box_plot, 
                  x = 'gene_rep', 
                  y = 'expr_dat',
                  fill = 'group_rep', 
                  xlab = FALSE,
                  ylab = 'Expression',
                  width = 0.6,#numeric value between 0 and 1 specifying box width.
                  bxp.errorbar = T, 
                  bxp.errorbar.width = 0.2, 
                  outlier.shape = NA, 
                  palette = 'npg' 
  ) + ylim(6.5, 11)
  p1 <- p1 + scale_x_discrete(
    breaks = c('PDCD1','CD274','CTLA4'),
    label = c('PD1','PD-L1','CTLA4'))
  p <- p1 + theme(
    plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 7.5, angle = 45,
                                 hjust = 0.5, vjust = 0.5),
    axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
    axis.title.x  = element_blank(),
    axis.title.y  = element_text(color = 'black', size = 10, angle = 90),
    legend.title  = element_blank(),
    legend.text   = element_text(color = 'black', size   = 7.5),
    axis.line = element_blank(), 
    panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA)
  )
  print(p)
  ggsave(filename = file.path(Figure_dir,'boxplot_riskScore_ICIgenes_expr_GSE65904_6-9.3cm.pdf'),
         plot = p, width = 6, height = 9.3, units = 'cm', dpi = 600)
  
  #Analyze the statistical significance
  if(T){
    p1 <- ggboxplot(box_plot, 
                    x = 'group_rep', 
                    y = 'expr_dat', 
                    fill = 'group_rep', 
                    xlab = FALSE,
                    ylab = 'Expression',
                    width = 0.6,#numeric value between 0 and 1 specifying box width.
                    bxp.errorbar = T, 
                    bxp.errorbar.width = 0.2, 
                    outlier.shape = NA, 
                    palette = 'npg' 
    ) 
    #Add significance markers
    p <-  p1 + geom_signif(comparisons = list(c("High_risk", "Low_risk")),
                           test = "t.test",
                           vjust=0.5,#move the text up or down relative to the bracket
                           ) + 
      facet_wrap(~gene_rep, nrow = 1,scales = 'fixed')
    print(p)
  }
}











