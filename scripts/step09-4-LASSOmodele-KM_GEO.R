#Survival comparison of risk score-high and -low group based on GEO dataset

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(survival)

#load exprSet and meta
load(file = file.path(Rdata_dir,paste(GSEName,'RData',sep=".")))
phe = na.omit(meta[,c(2,3)])
colnames(phe) = c('status_DSS','time_DSS')
phe = phe[(phe$time_DSS>=30),] 
phe$time_DSS = phe$time_DSS/365
dim(phe) 
exprSet = exprSet[,rownames(phe)]
dim(exprSet) 

#group patients based on RS
if(T){
  load(file = file.path(Rdata_dir,'step50_LASSO_modeleGenes_coef.RData'))
  LASSOgenes = rownames(coef_out3)
  exprSet_t = t(exprSet)
  expr_LASSOgenes = exprSet_t[,LASSOgenes]
  phe = cbind(phe,expr_LASSOgenes)
  phe$riskScore = coef_out3[1,1]*phe[,3] + coef_out3[2,1]*phe[,4] +
    coef_out3[3,1]*phe[,5] + coef_out3[4,1]*phe[,6] +
    coef_out3[5,1]*phe[,7] + coef_out3[6,1]*phe[,8] + coef_out3[7,1]*phe[,9]
  phe$RS_group = ifelse(phe$riskScore>median(phe$riskScore),'high','low')
  phe = phe[order(phe$riskScore,decreasing = FALSE),]
}
phe_sur = phe[,c('time_DSS','status_DSS','RS_group')]

#KM curve 
if(T){
  sfit <- survfit(Surv(time_DSS, status_DSS)~RS_group, data=phe_sur)
  sfit 
  summary(sfit)
  library(survminer)
  survplot = ggsurvplot(sfit,
                        title = 'Validation cohort', 
                        size = 0.4,
                        pval =TRUE,
                        pval.method = F, 
                        pval.size = 2.5, #text size of p
                        pval.coord = c(13,0.6), #position of p
                        xlab ="Survival time (Years)",
                        legend = c(0.8,0.94), 
                        legend.title = "", 
                        legend.labs = c("High risk", "Low risk"), # Change legend labels 
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
                          legend.title  = element_text(color = 'black', size  = 6),
                          legend.text   = element_text(color = 'black', size   = 6),
                          axis.ticks.x = element_line(size = 0.5),
                          axis.ticks.y = element_line(size = 0.5),
                          axis.line.y = element_line(color = 'black', linetype = 'solid',size = 0.4), 
                          axis.line.x = element_line (color = 'black',linetype = 'solid',size = 0.4), 
                        )
  )
  print(survplot)
  save(survplot, file = file.path(Rdata_dir,'step50_LASSOmodele_survplot_GEO.Rdata'))
  ggsave(filename = file.path(Figure_dir,'step50_LASSOmodele_KM_GEO_4.4-4.8cm.pdf'),
         plot = last_plot(), width = 4.4, height = 4.8, units = 'cm', dpi = 600)
  
  #Test for statistical significance with function survdiff（log-rank test）
  survdiff=survdiff(Surv(time_DSS, status_DSS)~RS_group,data=phe_sur)
  p.val = 1 - pchisq(survdiff$chisq, length(survdiff$n) - 1)
  p.val 
}



