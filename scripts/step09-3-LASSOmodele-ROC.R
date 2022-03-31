#Plot time-dependent ROC curve with survivalROC

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(purrr)
library(survivalROC)
library(ggplot2)

load(file = file.path(Rdata_dir,'step50_LASSO_KM_input.RData'))
phe_roc = phe[,c('time','status','riskScore')]
#3 columns are needed to plot ROC curve：syrvival time（year），survival status，riskScore

#1.plot ROC curves
if(T){
  if(T){
    cutoff = 1
    roc_1 = survivalROC(Stime=phe_roc$time, 
                        status=phe_roc$status, #1==dead，0==alive
                        marker = phe_roc$riskScore,
                        predict.time =cutoff, #predict.time：time point of the ROC curve
                        method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_1)

    cutoff = 3
    roc_3 = survivalROC(Stime=phe_roc$time, 
                        status=phe_roc$status, 
                        marker = phe_roc$riskScore,
                        predict.time =cutoff, #predict.time：time point of the ROC curve
                        method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_3)
    
    cutoff = 5
    roc_5 = survivalROC(Stime=phe_roc$time, 
                        status=phe_roc$status, 
                        marker = phe_roc$riskScore,
                        predict.time =cutoff, #predict.time：time point of the ROC curve
                        method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_5)
    
    cutoff = 10
    roc_10 = survivalROC(Stime=phe_roc$time, 
                        status=phe_roc$status, 
                        marker = phe_roc$riskScore,
                        predict.time =cutoff, #predict.time：time point of the ROC curve
                        method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_10)
    
    cutoff = 15
    roc_15 = survivalROC(Stime=phe_roc$time, 
                         status=phe_roc$status, 
                         marker = phe_roc$riskScore,
                         predict.time =cutoff, #predict.time：time point of the ROC curve
                         method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_15)
    
    roc1 = data.frame(roc_1$FP,roc_1$TP)
    colnames(roc1) = c('FP','TP')
    roc1$group = rep('1',nrow(roc1))
    roc3 = data.frame(roc_3$FP,roc_3$TP)
    colnames(roc3) = c('FP','TP')
    roc3$group = rep('3',nrow(roc3))
    roc5 = data.frame(roc_5$FP,roc_5$TP)
    colnames(roc5) = c('FP','TP')
    roc5$group = rep('5',nrow(roc5))
    roc10 = data.frame(roc_10$FP,roc_10$TP)
    colnames(roc10) = c('FP','TP')
    roc10$group = rep('10',nrow(roc10))
    roc15 = data.frame(roc_15$FP,roc_15$TP)
    colnames(roc15) = c('FP','TP')
    roc15$group = rep('15',nrow(roc15))
  }
 
  ## Plot
  roc = rbind(roc1,roc3,roc5)
  roc$TP = ifelse(roc$TP>1,1,roc$TP)
  if(T){
    p <- ggplot(roc, aes(x=FP, y=TP, colour=group,group=group)) + 
      geom_line(size=0.46)+
      geom_abline(intercept = 0,slope = 1,lty=2,lwd=0.3,colour="black")+
      xlim(0,1) + ylim(0,1) +
      labs(title="Training cohort",
           x="False positive rate",
           y="True positive rate") +
      scale_colour_manual(values = c('red','blue','green'),
                          limits = c('1','3','5'),
                          labels = c(paste("1 year AUC: ",round(roc_1$AUC,3)),
                                     paste("3 year AUC: ",round(roc_3$AUC,3)),
                                     paste("5 year AUC: ",round(roc_5$AUC,3)))
                          ) 
    p_roc <- p + theme(
      panel.background = element_rect(fill = 'white'), 
      legend.background = element_rect(fill = 'transparent'),
      legend.key = element_rect(fill = 'transparent'),
      panel.grid=element_blank(),
      plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
      plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
      plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
      axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
      axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
      axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
      axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
      axis.ticks.x = element_line(size = 0.5),
      axis.ticks.y = element_line(size = 0.5),
      legend.position=c(0.63,0.17),
      legend.key.size = unit(6.5, "pt"),
      legend.title  = element_blank(),
      legend.text   = element_text(color = 'black', size   = 5.8),
      axis.line = element_blank(), 
      panel.border = element_rect(linetype = 'solid', size = 0.4,fill = NA) 
    )
    print(p_roc)
    save(p_roc, file = file.path(Rdata_dir,'step50_LASSOmodele_p_roc.Rdata'))
    ggsave(filename = file.path(Figure_dir,'step50_LASSOmodele_ROC_1-3-5year_4.5-4.8cm.pdf'),
           plot = last_plot(), width = 4.5, height = 4.8, units = 'cm', dpi = 600)
  }
}





