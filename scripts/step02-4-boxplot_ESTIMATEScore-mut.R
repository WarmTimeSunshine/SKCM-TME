#Association between ESTIMATE Score and mutation

#rm(list=ls())
options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
figure_dir='../figures/'
table_dir = '../tables/'

score_f=file.path(Rdata_dir,
                  paste(paste("TCGA",request_cancer,'tumor_estimateScore.RData',sep="-")))
load(file = score_f)

#Read samples of wild and mutated
if(T){
  case_all_1 = read.table(file = file.path(table_dir,"case_mut-WT_1.tsv"),
                          header=T, fill = TRUE)
  case_all_2 = read.table(file = file.path(table_dir,"case_mut-WT_2.tsv"),
                          header=T, fill = TRUE)
  case_all_3 = read.table(file = file.path(table_dir,"case_mut-WT_3.tsv"),
                          header=T, fill = TRUE)
  case_all_4 = read.table(file = file.path(table_dir,"case_mut-WT_4.tsv"),
                          header=T, fill = TRUE)
  case_all_5 = read.table(file = file.path(table_dir,"case_mut-WT_5.tsv"),
                          header=T, fill = TRUE)
  case_all = rbind(case_all_1,case_all_2,case_all_3,case_all_4,case_all_5)
  case_all = case_all[,c(1,2)]
}

if(T){
  case_mut_1 = read.table(file = file.path(table_dir,"case_mut_1.tsv"),
                          header=T, fill = TRUE)
  case_mut_2 = read.table(file = file.path(table_dir,"case_mut_2.tsv"),
                          header=T, fill = TRUE)
  case_mut_3 = read.table(file = file.path(table_dir,"case_mut_3.tsv"),
                          header=T, fill = TRUE)
  case_mut = rbind(case_mut_1,case_mut_2,case_mut_3)
  case_mut = case_mut[,c(1,2)]
  case_mut$BRAF_type = rep("Mutant",times=nrow(case_mut))
}

case_wild = as.data.frame(case_all[!(rownames(case_all) %in% rownames(case_mut)),])
case_wild$BRAF_type = rep("WT",times=nrow(case_wild))
case = rbind(case_mut,case_wild)
case$submitter_id = rownames(case)

#To get dataframe mut_score_patient including mutation information and ESTIMATEScore
if(T){
  estimate_score_new = data.frame(estimate_score)
  estimate_score_new$submitter_id = substr(estimate_score_new$ID,1,12)
  estimate_score_new = estimate_score_new[!duplicated(estimate_score_new$submitter_id),]
  mut_score_patient = merge(case, estimate_score_new,by="submitter_id")
  mut_score_patient = mut_score_patient[,c(1,4,6:8)]
  dim(mut_score_patient)
}

#Boxplot 
if(T){
  require(cowplot)
  require(ggplot2)
  require(ggsci)
  require(ggpubr)
  dat=mut_score_patient
  
  #Stromal score---mutation
  if(T){
    p <- ggboxplot(dat, 
                   x = 'BRAF_type', 
                   y = 'Stromal_score', 
                   fill = 'BRAF_type', 
                   xlab = 'BRAF status', 
                   ylab = 'Stromal score',
                   bxp.errorbar = T, #error bar
                   bxp.errorbar.width = 0.2, 
                   outlier.shape = NA, #outlier
                   order = c('Mutant','WT'), 
                   palette = 'npg', 
                   legend = "none" 
    )+ ylim(-2000, 2000)

        my_comparisons <- list(c('Mutant','WT'))
    p <- p + geom_signif(comparisons = my_comparisons,
                         test = "t.test",
                         y_position	= 1700, #numeric vector with the y positions of the brackets
                         vjust=-0.2,#move the text up or down relative to the bracket
                         textsize=2) 

        p_boxplot_StromalScore_BRAF_mut <- p + theme(
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
    print(p_boxplot_StromalScore_BRAF_mut)
    save(p_boxplot_StromalScore_BRAF_mut, file = file.path(Rdata_dir,'p_boxplot_StromalScore_BRAF_mut.Rdata'))
    ggsave(filename = file.path(figure_dir,paste('fig2',request_cancer,'boxplot_StromalScore-BRAF_mut.pdf',sep = '_')),
           plot = last_plot(), width = 4.4, height = 4, units = 'cm', dpi = 600)
  }
  
  #Immune score---mutation
  if(T){
    p <- ggboxplot(dat, 
                   x = 'BRAF_type', 
                   y = 'Immune_score', 
                   fill = 'BRAF_type', 
                   xlab = 'BRAF status', 
                   ylab = 'Immune score',
                   bxp.errorbar = T, 
                   bxp.errorbar.width = 0.2, 
                   outlier.shape = NA, 
                   order = c('Mutant','WT'), 
                   palette = 'npg', 
                   legend = "none" 
    )+ ylim(-2000, 4500)

        my_comparisons <- list(c('Mutant','WT'))
    p <- p + geom_signif(comparisons = my_comparisons,
                         test = "t.test",
                         y_position	= 4000, #numeric vector with the y positions of the brackets
                         vjust=-0.2,#move the text up or down relative to the bracket
                         textsize=2) 
    
    p_boxplot_ImmuneScore_BRAF_mut <- p + theme(
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
    print(p_boxplot_ImmuneScore_BRAF_mut)
    save(p_boxplot_ImmuneScore_BRAF_mut, file = file.path(Rdata_dir,'p_boxplot_ImmuneScore_BRAF_mut.Rdata'))
    ggsave(filename = file.path(figure_dir,paste('fig2',request_cancer,'boxplot_ImmuneScore-BRAF_mut.pdf',sep = '_')),
           plot = last_plot(), width = 4.4, height = 4, units = 'cm', dpi = 600)
  }
  
}




