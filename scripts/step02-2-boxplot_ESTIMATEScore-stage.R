#Boxplot of association between ESTIMATE Score and stage

#rm(list=ls())
options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
figure_dir='../figures/'
table_dir = '../tables/'
data_type="Gene Expression Quantification" 
workflow_type="HTSeq - Counts"
#Load meta
exprSym_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,'exprSym_meta_groupT-N.RData',sep="-")))
load(file = exprSym_f)
#Read estimate_score
estimate_score = read.table(file = file.path(table_dir, paste("ESTIMATE",request_cancer,'fromWeb.txt',sep="_")),
                      header=TRUE)

#Step 1:Merge meta and estimate_score
# To get dataframe stage_score_patient including stage and ESTIMATEScore
if(T){
  id_stage=meta[,c('submitter_id','ajcc_pathologic_stage','tumor_stage')]
  id_stage=na.omit(id_stage)
  id_stage=id_stage[!(id_stage$ajcc_pathologic_stage %in% c("Not Reported",'Stage 0')),]
  dim(id_stage)
  library(stringr)
  id_stage$tumor_stage=str_split(id_stage$tumor_stage,' ',simplify = TRUE)[,2]
  id_stage$tumor_stage = toupper(id_stage$tumor_stage)
  id_stage$tumor_stage = gsub(pattern = "A", replacement = "", x = id_stage$tumor_stage)
  id_stage$tumor_stage = gsub(pattern = "B", replacement = "", x = id_stage$tumor_stage)
  id_stage$tumor_stage = gsub(pattern = "C", replacement = "", x = id_stage$tumor_stage)
  table(id_stage$tumor_stage)
  
  estimate_score = data.frame(estimate_score)
  estimate_score$ID = substr(estimate_score$ID,1,12)
  stage_score_patient = merge(id_stage, estimate_score,
                              by.x="submitter_id", by.y='ID',all=F)
  dim(stage_score_patient)
  colnames(stage_score_patient)=c('submitter_id','Pathologic_stage','Tumor_stage','Stromal_score',
                                  'Immune_score','ESTIMATE_score')
}

#Step 2: boxplot
if(T){
  require(cowplot)
  require(ggplot2)
  require(ggsci)
  require(ggpubr)
  dat=stage_score_patient
  
  #Boxplot of Stromal score and stage----
  #Determine differences between groups
  bartlett.test(Stromal_score~Tumor_stage,data = dat)
  shapiro.test(dat$Stromal_score)
  dat_anova<-aov(Stromal_score~Tumor_stage,data = dat) #ANOVA
  summary(dat_anova) 
  TukeyHSD(dat_anova)#post hoc test

  #Boxplot
  p <- ggboxplot(dat, 
                 x = 'Tumor_stage', 
                 y = 'Stromal_score', 
                 fill = 'Tumor_stage', 
                 xlab = 'Stage',
                 ylab = 'Stromal score',
                 bxp.errorbar = T, 
                 bxp.errorbar.width = 0.2, 
                 outlier.shape = NA, 
                 order = c('I','II','III','IV'), 
                 palette = 'npg', 
                 #add = 'point', 
                 legend = "none" 
                 )+ ylim(-2000, 2200)
  p <- p + theme(
    plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
    axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
    legend.title  = element_text(color = 'black', size  = 7.5),
    legend.text   = element_text(color = 'black', size   = 7.5),
    axis.line = element_blank(), 
    panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
  )
  print(p)
  dev.off()
  ggsave(filename = file.path(figure_dir,paste('fig2',request_cancer,'boxplot_StromalScore-Stage.pdf',sep = '_')),
    plot = p, width = 4.4, height = 4, units = 'cm', dpi = 300 )
  
  #Boxplot of Immune score and stage----
  bartlett.test(Immune_score~Tumor_stage,data = dat) 
  shapiro.test(dat$Immune_score)
  dat_anova<-aov(Immune_score~Tumor_stage,data = dat) 
  summary(dat_anova) 
  TukeyHSD(dat_anova)

  p <- ggboxplot(dat, 
                 x = 'Tumor_stage', 
                 y = 'Immune_score', 
                 fill = 'Tumor_stage', 
                 xlab = 'Stage',
                 ylab = 'Immune score',
                 bxp.errorbar = T, 
                 bxp.errorbar.width = 0.2, 
                 outlier.shape = NA, 
                 order = c('I','II','III','IV'),
                 palette = 'npg', 
                 legend = "none" 
  ) + ylim(-2000, 4800)

  my_comparisons <- list( c("I", "II"), c("II", "III"))
  p <-  p + stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                               vjust = 0.5,
                               symnum.args = list(cutpoints = c(0, 0.05, 1), symbols = c("*", "ns")) 
  )
  p <- p + theme(
    plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
    axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
    legend.title  = element_text(color = 'black', size  = 7.5),
    legend.text   = element_text(color = 'black', size   = 7.5),
    axis.line = element_blank(), 
    panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
  )
  print(p)
  dev.off()
  ggsave(filename = file.path(figure_dir,paste('fig2',request_cancer,'boxplot_ImmuneScore-Stage.pdf',sep = '_')),
         plot = p, width = 4.4, height = 4, units = 'cm', dpi = 600 )
}





