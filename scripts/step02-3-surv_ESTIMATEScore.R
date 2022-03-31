#KM analysis between immune score-high and -low patients

library(survival)
library(survminer)
#rm(list=ls())
options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
figure_dir='../figures/'
table_dir = '../tables/'
data_type="Gene Expression Quantification" 
workflow_type="HTSeq - Counts"
#load meta
exprSym_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,'exprSym_meta_groupT-N.RData',sep="-")))
load(file = exprSym_f)
estimate_score = read.table(file = file.path(table_dir, paste("ESTIMATE",request_cancer,'fromWeb.txt',sep="_")),
                            header=TRUE)

#To get dataframe surv_scoreGroup_patient including survival information and immune score
if(T){
  clin_info = meta[,c('submitter_id','vital_status','days_to_death','days_to_last_follow_up')]
  clin_info = as.data.frame(clin_info)
  clin_info$time = ifelse(clin_info$vital_status=="Alive",
                          clin_info$days_to_last_follow_up, clin_info$days_to_death)
  clin_info = clin_info[!is.na(clin_info$time),]
  clin_info = clin_info[clin_info$time>=30,]
  clin_info$time = clin_info$time/365
  clin_info$status=ifelse(clin_info$vital_status=="Dead",1,0)

  estimate_score = data.frame(estimate_score)
  estimate_score$ID = substr(estimate_score$ID,1,12)
  dup_samples = estimate_score$ID[duplicated(estimate_score$ID)]
  estimate_score=estimate_score[!(estimate_score$ID %in% dup_samples),]
  
  surv_score_patient = merge(clin_info, estimate_score,
                             by.x="submitter_id", by.y='ID')
  dim(surv_score_patient)
  surv_scoreGroup_patient = surv_score_patient[,c('submitter_id','vital_status','time',
                                                  'status','Stromal_score','Immune_score','ESTIMATE_score')]
  surv_scoreGroup_patient = na.omit(surv_scoreGroup_patient)
  dim(surv_scoreGroup_patient)
  
  #Divide patients based on median 
  surv_scoreGroup_patient$StromalScore_group = ifelse(surv_scoreGroup_patient$Stromal_score
                                                 >median(surv_scoreGroup_patient$Stromal_score),
                                                 "high","low")
  surv_scoreGroup_patient$ImmuneScore_group = ifelse(surv_scoreGroup_patient$Immune_score
                                                >median(surv_scoreGroup_patient$Immune_score),
                                                "high","low")
}

#KM curve based on Stromal score 
if(T){
  sfit <- survfit(Surv(time, status)~StromalScore_group, data=surv_scoreGroup_patient)
  sfit 
  summary(sfit)
  survplot = ggsurvplot(sfit,
                        title = 'Stromal score',
                        size = 0.4,
                        pval =TRUE,
                        pval.method = F, 
                        pval.size = 2.5, #text size of p
                        pval.coord = c(22,0.5), #position of p
                        xlab ="Survival time (Years)",
                        legend = c(0.8,0.9), 
                        legend.title = "", 
                        legend.labs = c("High", "Low"), # Change legend labels 
                        ggtheme = theme(
                          panel.background = element_rect(fill = 'white'), 
                          legend.key = element_rect(fill = 'white'),
                          plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
                          plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
                          plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
                          axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
                          axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
                          axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
                          axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
                          axis.ticks.x = element_line(size = 0.5),
                          axis.ticks.y = element_line(size = 0.5),
                          legend.title  = element_text(color = 'black', size  = 6),
                          legend.text   = element_text(color = 'black', size   = 6),
                          axis.line.y = element_line(color = 'black', linetype = 'solid',size = 0.4), 
                          axis.line.x = element_line (color = 'black',linetype = 'solid',size = 0.4), 
                        )
  ) 
  print(survplot)
  ggsave(filename = file.path(figure_dir,paste('fig2',request_cancer,'surv_StromalScore.pdf',sep = '_')),
         plot = last_plot(), width = 4.4, height = 4.8, units = 'cm', dpi = 600)
  
  #Test for statistical significance（log-rank test）
  survdiff=survdiff(Surv(time, status)~StromalScore_group,data=surv_scoreGroup_patient)
  p.val = 1 - pchisq(survdiff$chisq, length(survdiff$n) - 1)
  p.val 
}

#KM curve based on Immune score 
if(T){
  sfit <- survfit(Surv(time, status)~ImmuneScore_group, data=surv_scoreGroup_patient)
  sfit 
  summary(sfit)
  survplot = ggsurvplot(sfit,
                        title = 'Immune score',
                        size = 0.4,
                        pval =TRUE,
                        pval.method = F, 
                        pval.size = 2.5, #text size of p
                        pval.coord = c(17,0.6), #position of p
                        xlab ="Survival time (Years)",
                        legend = c(0.8,0.94), 
                        legend.title = "", 
                        legend.labs = c("High", "Low"), # Change legend labels 
                        ggtheme = theme(
                          panel.background = element_rect(fill = 'white'), 
                          legend.key = element_rect(fill = 'white'),
                          plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
                          plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
                          plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
                          axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
                          axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
                          axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
                          axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
                          axis.ticks.x = element_line(size = 0.5),
                          axis.ticks.y = element_line(size = 0.5),
                          legend.title  = element_text(color = 'black', size  = 6),
                          legend.text   = element_text(color = 'black', size   = 6),
                          axis.line.y = element_line(color = 'black', linetype = 'solid',size = 0.4), 
                          axis.line.x = element_line (color = 'black',linetype = 'solid',size = 0.4), 
                        )
  ) 
  print(survplot)
  ggsave(filename = file.path(figure_dir,paste('fig2',request_cancer,'surv_ImmuneScore.pdf',sep = '_')),
         plot = last_plot(), width = 4.4, height = 4.8, units = 'cm', dpi = 600)
  
  survdiff=survdiff(Surv(time, status)~ImmuneScore_group,data=surv_scoreGroup_patient)
  p.val = 1 - pchisq(survdiff$chisq, length(survdiff$n) - 1)
  p.val 
}

#Correlation curve：Stromal score VS immune score
if(T){
  cor(estimate_score$Stromal_score,estimate_score$Immune_score,method = "pearson")
  library(ggpubr)
  p <- ggscatter(estimate_score,
            x = "Stromal_score", 
            y = "Immune_score",
            color = 'black', #color of points
            size =0.5,  #size of points
            add = "reg.line",
            add.params = list(color = "cyan4",size = 0.9),
            conf.int = F,
            cor.coef = TRUE, 
            cor.method = "pearson",
            cor.coef.size = 2.5,
            xlab = "Stromal score", 
            ylab = "Immune score",
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
  )+ylim(-2000,4500)
  print(p)
  dev.off()
  ggsave(filename = file.path(figure_dir,paste('fig2',request_cancer,'cor_Stromal-ImmuneScore.pdf',sep = '_')),
         plot = last_plot(), width = 5.2, height = 4.8, units = 'cm', dpi = 600)
}




