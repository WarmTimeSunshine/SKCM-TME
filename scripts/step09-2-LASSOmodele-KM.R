#Survival comparison of risk score-high and -low group

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(survival)

load(file = file.path(Rdata_dir,'step50_LASSO_KM_input.RData'))
phe_sur = phe[,c('time','status','RS_group')]

#KM curve 
sfit <- survfit(Surv(time, status)~RS_group, data=phe_sur)
sfit 
summary(sfit)
library(survminer)
survplot = ggsurvplot(sfit,
                      title = 'Training cohort', 
                      size = 0.4,
                      pval =TRUE,
                      pval.method = F, 
                      pval.size = 2.5, #text size of p
                      pval.coord = c(19,0.58), #position of p
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
save(survplot, file = file.path(Rdata_dir,'step50_LASSOmodele_survplot.Rdata'))
ggsave(filename = file.path(Figure_dir,'step50_LASSOmodele_KM_4.4-4.8cm.pdf'),
       plot = last_plot(), width = 4.4, height = 4.8, units = 'cm', dpi = 600)

#Test for statistical significance with function survdiff（log-rank test）
survdiff=survdiff(Surv(time, status)~RS_group,data=phe_sur)
p.val = 1 - pchisq(survdiff$chisq, length(survdiff$n) - 1)
p.val


