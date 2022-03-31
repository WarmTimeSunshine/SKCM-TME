#To analyze the association between RS and immune infiltration（from ImmuCellAI website）

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(ggplot2)
library(ggpubr)

#step 01:Process expression data and group information to input into ImmuneCellAI website
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
  
  #group patients
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
                                    'high_risk','low_risk')
  samGroup = (t(expr_lassoGenes))['RS_group',]
  samRS = (t(expr_lassoGenes))['riskScore',]
  identical(colnames(exprDat),names(samGroup))
  exprDat = rbind(samGroup,exprDat)
  rownames(exprDat)[1] = 'group'
  expr_to_write = rbind(colnames(exprDat),exprDat)
  expr_to_write = cbind(rownames(expr_to_write),expr_to_write)
  expr_to_write[1,1] = 'Symbol'
  write.table(expr_to_write, file=file.path(Table_dir,'ImmuCellAI_SKCM_expr.txt'),
              sep = '\t', 
              quote = F,
              row.names=F, col.names=F)
}

#step 02:Correlation of riskScore and immune infiltration in TCGA dataset
if(T){
  immune_infil = read.table(file = file.path(Table_dir,'ImmuCellAI_SKCM_result_immune_infilatration.txt'),
                            header = T)
  #To get the matrix dat_plot for boxplot
  if(T){
    colnames(immune_infil) = c('Response','CD4_naive','CD8_naive','Tc',
                               'Tex','Tr1','nTreg','iTreg',
                               'Th1','Th2','Th17','Tfh',
                               'Tcm','Tem','NKT','MAIT',
                               'DC','B_cell','Monocyte','Macrophage',
                               'NK','Neutrophil','Tgd','CD4+T','CD8+T',
                               'InfiltrationScore')
    immune_plot = immune_infil[,2:25]
    sample_rep = rep(rownames(immune_plot),ncol(immune_plot))
    identical(names(samRS),rownames(immune_plot))
    RS_rep = rep(samRS,ncol(immune_plot))
    z = array(colnames(immune_plot),dim=c(24,1))
    cell_matrix = apply(z, 1, function(x) rep(x,nrow(immune_plot)))
    cell_rep = as.vector(cell_matrix)
    immune_dat = as.vector(as.matrix(immune_plot))
    
    dat_plot = as.data.frame(cbind(sample_rep,cell_rep,RS_rep,immune_dat))
    dat_plot$immune_dat = as.numeric(dat_plot$immune_dat)
    dat_plot$RS_rep = as.numeric(dat_plot$RS_rep)
    dat_plot$cell_rep = factor(dat_plot$cell_rep,levels = c('CD4+T','CD8+T','CD4_naive','CD8_naive','Tc',
                                                            'Tex','Tcm','Tem',
                                                            'Tr1','nTreg','iTreg',
                                                            'Th1','Th2','Th17','Tfh',
                                                            'MAIT','Tgd','NKT',
                                                            'NK','DC','B_cell','Macrophage','Monocyte',
                                                            'Neutrophil'))
  }
  
  library(ggpubr)
  p1 <- ggscatter(dat_plot,
                  x = "RS_rep", 
                  y = "immune_dat",
                  color = 'black', #color of points
                  size =0.2,  #size of points
                  add = "reg.line",
                  add.params = list(color = "blue",size=0.6),
                  conf.int = F,
                  cor.coef = TRUE, 
                  cor.method = "pearson",
                  cor.coef.size = 2,
                  xlab = "Risk Score", 
                  ylab = "Abundance",
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
  )+ylim(0.0,0.8)
  p <- p1+facet_wrap(.~cell_rep, 
                     scales = 'fixed',
                     nrow =4 
  )+theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text  = element_text(color = 'black',size = 7.5, angle = 0)
  )
  print(p)
  ggsave(filename = file.path(Figure_dir,'cor_riskScore-ImmuneInfil_17-12cm.pdf'),
         plot = last_plot(), width = 17, height = 12, 
         units = 'cm', dpi = 600)
}

#step 03:Boxplot of riskScore and immune infiltration in TCGA dataset
if(T){
  immune_infil = read.table(file = file.path(Table_dir,'ImmuCellAI_SKCM_result_immune_infilatration.txt'),
                            header = T)
  if(T){
    colnames(immune_infil) = c('Response','CD4_naive','CD8_naive','Tc',
                               'Tex','Tr1','nTreg','iTreg',
                               'Th1','Th2','Th17','Tfh',
                               'Tcm','Tem','NKT','MAIT',
                               'DC','B_cell','Monocyte','Macrophage',
                               'NK','Neutrophil','Tgd','CD4+T','CD8+T',
                               'InfiltrationScore')
    immune_plot = immune_infil[,2:25]
    sample_rep = rep(rownames(immune_plot),ncol(immune_plot))
    identical(names(samGroup),rownames(immune_plot))
    group_rep = rep(samGroup,ncol(immune_plot))
    z = array(colnames(immune_plot),dim=c(24,1))
    cell_matrix = apply(z, 1, function(x) rep(x,nrow(immune_plot)))
    cell_rep = as.vector(cell_matrix)
    immune_dat = as.vector(as.matrix(immune_plot))
    
    dat_plot = as.data.frame(cbind(sample_rep,cell_rep,group_rep,immune_dat))
    dat_plot$immune_dat = as.numeric(dat_plot$immune_dat)
  }
  
  #1.boxplot
  if(T){
    #data procesing
    if(T){
      dat = dat_plot
      dat$cell_rep = factor(dat$cell_rep,levels = c('CD4+T','CD8+T','CD4_naive','CD8_naive','Tc',
                                                    'Tex','Tcm','Tem',
                                                    'Tr1','nTreg','iTreg',
                                                    'Th1','Th2','Th17','Tfh',
                                                    'MAIT','Tgd','NKT',
                                                    'NK','DC','B_cell','Macrophage','Monocyte',
                                                    'Neutrophil'))
      cell_type = c('CD4+T','CD8+T','CD4_naive','CD8_naive','Tc',
                    'Tex','Tcm','Tem',
                    'Tr1','nTreg','iTreg',
                    'Th1','Th2','Th17','Tfh',
                    'MAIT','Tgd','NKT',
                    'NK','DC','B_cell','Macrophage','Monocyte',
                    'Neutrophil')
      dat$group_rep = sapply(strsplit(dat$group_rep,'_'),"[[",1)
      dat$group_rep = capitalize(dat$group_rep)
    }
    #plot in sequence
    if(T){
      p_all = list()
      i=1
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p <- ggboxplot(dt, 
                       x = 'group_rep', 
                       y = 'immune_dat', 
                       fill = 'group_rep',
                       xlab = NULL,
                       ylab = NULL,
                       title = cell_type[i],
                       width = 0.6,#numeric value between 0 and 1 specifying box width.
                       bxp.errorbar = T, 
                       bxp.errorbar.width = 0.2, 
                       outlier.shape = NA, 
                       palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.03)
        p = p + geom_signif(comparisons = list(c("High", "Low")),
                            test = "t.test",
                            y_position	= 0.34, #numeric vector with the y positions of the brackets
                            vjust=-0.5,#move the text up or down relative to the bracket
                            textsize=2.4)
        p <- p + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p)
        p_all[[i]] = p
      }
      i=2
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p2 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.1)
        p2 = p2 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.45, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p2 <- p2 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p2)
        p_all[[i]] = p2
      }
      i=3
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p3 <- ggboxplot(dt, 
                        x = 'group_rep',
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, #error bar
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.07)
        p3 = p3 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position=0.24, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p3 <- p3 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p3)
        p_all[[i]] = p3
      }
      i=4
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat',
                        fill = 'group_rep',
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.03)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.24, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none",
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=5
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.14)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.62, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=6
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.1)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.35, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=7
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2,
                        outlier.shape = NA, 
                        palette = 'npg'
        ) + ylim(0, max(dt$immune_dat)+0.07)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.24, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=8
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.03)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none",
          axis.line.y = element_line(color = 'black', linetype = 'solid'),
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=9
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA,
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.07)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=10
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg'
        ) + ylim(0, max(dt$immune_dat)+0.00)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.15, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=11
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep',
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.1)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.37, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=12
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep',
                        y = 'immune_dat',
                        fill = 'group_rep',
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2,
                        outlier.shape = NA,
                        palette = 'npg'
        ) + ylim(0, max(dt$immune_dat)+0.1)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.38, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'),
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=13
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep',
                        y = 'immune_dat', 
                        fill = 'group_rep',
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg'
        ) + ylim(0, max(dt$immune_dat)+0.1)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.5, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none",
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=14
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.025)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.37, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none",
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=15
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T,
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.146)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.58, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none",
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=16
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA,
                        palette = 'npg'
        ) + ylim(0, max(dt$immune_dat)+0.11)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.45, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'),
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=17
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T,
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.1)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.37, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=18
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.1)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.33, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'),
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=19
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep',
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T,
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.125)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.52, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'),
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=20
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.152)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.6, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=21
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)-0.1)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.5, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=22
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat',
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T,
                        bxp.errorbar.width = 0.2,
                        outlier.shape = NA, 
                        palette = 'npg'
        ) + ylim(0, max(dt$immune_dat)+0.2)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.75, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'),
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=23
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat',
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T,
                        bxp.errorbar.width = 0.2,
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.085)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              #y_position	= , #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
      i=24
      if(T){
        dt = dat[dat$cell_rep==cell_type[i],]
        p1 <- ggboxplot(dt, 
                        x = 'group_rep', 
                        y = 'immune_dat', 
                        fill = 'group_rep', 
                        xlab = NULL,
                        ylab = NULL,
                        title = cell_type[i],
                        width = 0.6,#numeric value between 0 and 1 specifying box width.
                        bxp.errorbar = T, 
                        bxp.errorbar.width = 0.2, 
                        outlier.shape = NA, 
                        palette = 'npg' 
        ) + ylim(0, max(dt$immune_dat)+0.12)
        p1 = p1 + geom_signif(comparisons = list(c("High", "Low")),
                              test = "t.test",
                              y_position	= 0.46, #numeric vector with the y positions of the brackets
                              vjust=-0.5,#move the text up or down relative to the bracket
                              textsize=2.4)
        p1 <- p1 + theme(
          plot.title    = element_text(color = 'black', size   = 7.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 7.5),
          legend.position = "none", 
          axis.line.y = element_line(color = 'black', linetype = 'solid'), 
          axis.line.x = element_line (color = 'black',linetype = 'solid'), 
        )
        print(p1)
        p_all[[i]] = p1
      }
    }
    #Combine images with R package patchwork
    if(T){
      library(patchwork)
      wrap_plots(p_all,nrow = 4,guides = "collect")
      ggsave(filename = file.path(Figure_dir,'boxplot_ImmuneInfil_all_TCGA.pdf'),
             plot = last_plot(), width = 17, 
             height = 15, 
             units = 'cm', dpi = 600)
    }
  }
  
  #2.BoxPlot of immune cells enriched in low-risk group
  if(T){
    if(T){
      dat = dat_plot[dat_plot$cell_rep %in% c('CD8+T','Tc','Tex','Tr1','iTreg',
                                              'Th1','Th2','Tfh','MAIT',
                                              'NK','DC','B_cell','Macrophage'),]
      dat$cell_rep = factor(dat$cell_rep,levels = c('CD8+T','Tc','Tex','Tr1','iTreg',
                                                    'Th1','Th2','Tfh','MAIT',
                                                    'NK','DC','B_cell','Macrophage'))
      library(Hmisc)
      dat$group_rep = capitalize(dat$group_rep)
    }
    
    if(T){
      p1 <- ggboxplot(dat, 
                      x = 'cell_rep', 
                      y = 'immune_dat', 
                      fill = 'group_rep',
                      xlab = 'Immune Cell',
                      ylab = 'Abundance',
                      width = 0.6,#numeric value between 0 and 1 specifying box width.
                      bxp.errorbar = T, 
                      bxp.errorbar.width = 0.2,
                      outlier.shape = NA, 
                      palette = 'npg' 
      ) + ylim(0, 0.75)
      p <- p1 + theme(
        plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 7.5, angle = 45,
                                     hjust = 0.5, vjust = 0.5),
        axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 10, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 10, angle = 90),
        legend.title  = element_blank(),
        legend.text   = element_text(color = 'black', size   = 7.5),
        legend.position ="top",
        axis.line = element_blank(), 
        panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
      )
      print(p)
      ggsave(filename = file.path(Figure_dir,'boxplot_ImmuneInfil_lowRS_TCGA_17-9.9cm.pdf'),
             plot = p, width = 17, height = 9.9, units = 'cm', dpi = 600)
      dev.off()
    }
    
    if(T){
      p1 <- ggboxplot(dat, 
                      x = 'group_rep', 
                      y = 'immune_dat',
                      width = 0.6,#numeric value between 0 and 1 specifying box width.
                      bxp.errorbar = T, 
                      bxp.errorbar.width = 0.2, 
                      outlier.shape = NA, 
                      order = c('High_risk','Low_risk'), 
                      palette = 'npg' 
      ) + ylim(0, 0.8)
      # #Add significance markers
      p <-  p1 + geom_signif(comparisons = list(c("High_risk", "Low_risk")),
                             test = "t.test",
                             vjust=0.5,#move the text up or down relative to the bracket
                             map_signif_level=c("*"=0.05)) + 
        facet_wrap(~cell_rep, nrow = 1,scales = 'fixed')
      p <- p + theme(
        plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 7.5, angle = 45,
                                     hjust = 0.5, vjust = 0.5),
        axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 10, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 10, angle = 90),
        legend.title  = element_blank(),
        legend.text   = element_text(color = 'black', size   = 7.5),
        legend.position ="top",
        axis.line.y = element_line(color = 'black', linetype = 'solid'), 
        axis.line.x = element_line (color = 'black',linetype = 'solid'), 
      )
      print(p)
      ggsave(filename = file.path(Figure_dir,'boxplot_ImmuneInfil_lowRS_sigSign.pdf'),
             plot = p, width = 17, height = 9.9, units = 'cm', dpi = 600)
      dev.off()
    }
  } 
  
  #3.Boxplot of immune cells enriched in high-risk group
  if(T){
    if(T){
      dat = dat_plot[dat_plot$cell_rep %in% c('CD4_naive','CD8_naive','Tcm',
                                              'Th17','NKT','Neutrophil'),]
      dat$cell_rep = factor(dat$cell_rep,levels = c('CD4_naive','CD8_naive','Tcm',
                                                    'Th17','NKT','Neutrophil'))
      library(Hmisc)
      dat$group_rep = capitalize(dat$group_rep)
    }
    
    if(T){
      p1 <- ggboxplot(dat, 
                      x = 'cell_rep', 
                      y = 'immune_dat', 
                      fill = 'group_rep', 
                      xlab = 'Immune Cell',
                      ylab = 'Abundance',
                      width = 0.6,#numeric value between 0 and 1 specifying box width.
                      bxp.errorbar = T, 
                      bxp.errorbar.width = 0.2,
                      outlier.shape = NA, 
                      palette = 'npg' 
      ) + ylim(0, 0.47)
      p <- p1 + theme(
        plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 7.5, angle = 45,
                                     hjust = 0.5, vjust = 0.5),
        axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 10, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 10, angle = 90),
        legend.title  = element_blank(),
        legend.text   = element_text(color = 'black', size   = 7.5),
        legend.position ="top",
        axis.line = element_blank(), 
        panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
      )
      print(p)
      ggsave(filename = file.path(Figure_dir,'boxplot_ImmuneInfil_HighRS_TCGA_9-9.9.pdf'),
             plot = p, width = 9, height = 9.9, units = 'cm', dpi = 600)
      dev.off()
    }
    
    if(T){
      p1 <- ggboxplot(dat, 
                      x = 'group_rep', 
                      y = 'immune_dat', 
                      width = 0.6,#numeric value between 0 and 1 specifying box width.
                      bxp.errorbar = T,
                      bxp.errorbar.width = 0.2,
                      outlier.shape = NA,
                      order = c('High_risk','Low_risk'), 
                      palette = 'npg' 
      ) + ylim(0, 0.6)
      # #Add significance markers
      p <-  p1 + stat_compare_means(comparisons = list(c("High_risk", "Low_risk")),label = "p.signif",
                                    symnum.args = list(cutpoints = c(0, 0.05, 1), symbols = c("*", "ns")) 
      ) + facet_wrap(~cell_rep, nrow = 1,scales = 'fixed')
      p <- p + theme(
        plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 7.5, angle = 45,
                                     hjust = 0.5, vjust = 0.5),
        axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 10, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 10, angle = 90),
        legend.title  = element_blank(),
        legend.text   = element_text(color = 'black', size   = 7.5),
        legend.position ="top",
        axis.line.y = element_line(color = 'black', linetype = 'solid'), 
        axis.line.x = element_line (color = 'black',linetype = 'solid'), 
      )
      print(p)
      ggsave(filename = file.path(Figure_dir,'boxplot_ImmuneInfil_HighRS_sigSign.pdf'),
             plot = p, width = 9, height = 9.9, units = 'cm', dpi = 600)
      dev.off()
    }
  }
}

#step 04:Correlation of RS with specific genes expression（PD1,PDL1,CTLA4）
if(T){
  #1.collect required matrix dat_plot
  if(T){
    x = exprSym_FPKM[,group_list_TN=="Tumor"]
    expr_log = log2(x+1)
    x=t(expr_log[c('CD274','PDCD1','CTLA4','MAGEA3'),])
    identical(rownames(x),rownames(expr_lassoGenes))
    expr = cbind(expr_lassoGenes, x)
    cor(expr$riskScore,expr[,c('CD274','PDCD1','CTLA4','MAGEA3')],method = "pearson")

    expr_plot = expr[,c(8,10,11,12)]
    sample_rep = rep(rownames(expr_plot),ncol(expr_plot)-1)
    RS_rep = rep(expr_plot$riskScore,ncol(expr_plot)-1)
    z = array(colnames(expr_plot[,-1]),dim=c(ncol(expr_plot[,-1]),1))
    gene_matrix = apply(z, 1, function(x) rep(x,nrow(expr_plot)))
    gene_rep = as.vector(gene_matrix)
    expr_dat = as.vector(as.matrix(expr_plot[,-1]))
    
    dat_plot = as.data.frame(cbind(sample_rep,RS_rep,gene_rep,expr_dat))
    dat_plot$expr_dat = as.numeric(dat_plot$expr_dat)
    dat_plot$RS_rep = as.numeric(dat_plot$RS_rep)
    dat_plot$gene_rep = factor(dat_plot$gene_rep,levels = c('CD274','PDCD1','CTLA4'))
  }
  
  #2.plot correlation curves in sequence
  if(T){
    gene_type = c('CD274','PDCD1','CTLA4')
    if(T){
      i=1
      dt = dat_plot[dat_plot$gene_rep==gene_type[1],]
      library(ggpubr)
      p <- ggscatter(dt,
                     x = "RS_rep", 
                     y = "expr_dat",
                     color = 'black', #color of points
                     size =0.4,  #size of points
                     add = "reg.line",
                     add.params = list(color = "blue",size=0.6),
                     conf.int = T,
                     cor.coef = TRUE,
                     cor.method = "pearson",
                     cor.coef.size = 2.5,
                     xlab = "Risk Score", 
                     ylab = paste0(gene_type[i],"  Expression"),
                     ggtheme =  theme(
                       panel.background = element_rect(fill = 'white'), 
                       legend.key = element_rect(fill = 'white'),
                       plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
                       plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
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
      )+ylim(0.0,5.5)
      print(p)
      ggsave(filename = file.path(Figure_dir,'cor_riskScore_expr_PD1_4.5-4.5cm.pdf'),
             plot = last_plot(), width = 4.5, 
             height = 4.5, 
             units = 'cm', dpi = 600)
    }
    if(T){
      i=2
      dt = dat_plot[dat_plot$gene_rep==gene_type[2],]
      library(ggpubr)
      p <- ggscatter(dt,
                     x = "RS_rep", 
                     y = "expr_dat",
                     color = 'black', #color of points
                     size =0.4,  #size of points
                     add = "reg.line",
                     add.params = list(color = "blue",size=0.6),
                     conf.int = T,
                     cor.coef = TRUE, 
                     cor.method = "pearson",
                     cor.coef.size = 2.5,
                     xlab = "Risk Score", 
                     ylab = paste0(gene_type[i],"  Expression"),
                     ggtheme =  theme(
                       panel.background = element_rect(fill = 'white'), 
                       legend.key = element_rect(fill = 'white'),
                       plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
                       plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
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
      )+ylim(0.0,6.5)
      print(p)
      ggsave(filename = file.path(Figure_dir,'cor_riskScore_expr_PDL1_4.5-4.5cm.pdf'),
             plot = last_plot(), width = 4.5, 
             height = 4.5, 
             units = 'cm', dpi = 600)
    }
    if(T){
      i=3
      dt = dat_plot[dat_plot$gene_rep==gene_type[3],]
      library(ggpubr)
      p <- ggscatter(dt,
                     x = "RS_rep", 
                     y = "expr_dat",
                     color = 'black', #color of points
                     size =0.4,  #size of points
                     add = "reg.line",
                     add.params = list(color = "blue",size=0.6),
                     conf.int = T,
                     cor.coef = TRUE,
                     cor.method = "pearson",
                     cor.coef.size = 2.5,
                     xlab = "Risk Score",
                     ylab = paste0(gene_type[i],"  Expression"),
                     ggtheme =  theme(
                       panel.background = element_rect(fill = 'white'), 
                       legend.key = element_rect(fill = 'white'),
                       plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
                       plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
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
      )+ylim(0.0,7)
      print(p)
      ggsave(filename = file.path(Figure_dir,'cor_riskScore_expr_CTLA4_4.5-4.5cm.pdf'),
             plot = last_plot(), width = 4.5, 
             height = 4.5, 
             units = 'cm', dpi = 600)
    }
  }
}














