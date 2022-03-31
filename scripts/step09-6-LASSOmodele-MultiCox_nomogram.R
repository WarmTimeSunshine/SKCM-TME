#Construct model based on multivariate Cox regression analysis, and plot forest map.
#Then construct nomogram

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library("survival")
library("survminer")
library('forestplot')

#1.load clinical data and riskScore
if(T){
  load(file = file.path(Rdata_dir,'step50_LASSO_KM_input.RData'))
  #phe includes survival time,survival status and riskScore
  
  workflow_type = "HTSeq - FPKM"
  exprSym_f=file.path(Rdata_dir,
                      paste(paste("TCGA",request_cancer,workflow_type,'exprSym_meta_groupT-N.RData',sep="-")))
  load(file = exprSym_f)
  dat_meta = meta[match(rownames(phe),meta$submitter_id),c('submitter_id','age_at_diagnosis','gender',
                                                           'ajcc_pathologic_stage','ajcc_pathologic_t',
                                                           'ajcc_pathologic_n','ajcc_pathologic_m')]
  identical(rownames(phe),dat_meta$submitter_id)
  dat_meta = cbind(phe[,c('time','status','riskScore')],dat_meta)
  dat_meta = na.omit(dat_meta[,-4])
}

#2.Data type conversion (classification variables such as gender, grade variables such as stage should be converted to numbers)  
if(T){
  dat=dat_meta
  dat$gender = ifelse(dat$gender=='male',1,2) #1 male，2 female
  
  colnames(dat)[4] = 'age_year'
  dat$age_year = dat$age_year/365 
  dat$age[dat$age_year<75] <- 1
  dat$age[dat$age_year>=75] <- 2
  
  #construct function to convert stage to number
  map_To_num <- function(x){
    a=10
    if(x %in% c('Stage 0','T0','N0','M0')) 
      a=1
    if(x %in% c('Stage I','Stage IA','Stage IB','Tis','T1','T1a','T1b',
                'N1','N1a','N1b','M1a','M1b','M1c')) 
      a=2
    if(x %in% c('Stage II','Stage IIA','Stage IIB','Stage IIC','T2','T2a','T2b',
                'N2','N2a','N2b','N2c')) 
      a=3
    if(x %in% c('Stage III','Stage IIIA','Stage IIIB','Stage IIIC','T3','T3a','T3b',
                'N3')) 
      a=4
    if(x %in% c('Stage IV','T4','T4a','T4b')) 
      a=5
    return(a)
  }
  dat_stage=apply(as.matrix(dat[,6:9]),c(1,2),map_To_num) 
  colnames(dat_stage) = c('stage','T','N','M')
  identical(rownames(dat),rownames(dat_stage))
  dat = cbind(dat,dat_stage)
}

#3.Univariate Cox regression
if(T){
  #Univariate Cox regression, result is saved as cox_results
  if(T){
    var_cox = dat[,c('age','gender','stage','T','N','M','riskScore')]
    sur_cox = dat[,c('time','status')]
    cox_results <-apply(var_cox, 2, function(var){
      sur_cox$varname = var
      m=coxph(Surv(time, status) ~ varname, data=sur_cox)
      beta <- coef(m) 
      se <- sqrt(diag(vcov(m)))
      HR <- exp(beta)
      HRse <- HR * se
      tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                         HR = HR, HRse = HRse,
                         HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                         HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                         HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
      return(tmp['varname',])
    })
    cox_results=t(cox_results)
  }
  
  dat_plot = as.data.frame(cox_results[,c('p','HR','HRCILL','HRCIUL')])
  dat_plot = cbind(rownames(cox_results),dat_plot)
  colnames(dat_plot)[1] = 'Varnames'
  library(Hmisc)
  dat_plot[,1] = capitalize(dat_plot[,1])
  
  #forestplot
  if(T){
    #1.construct data txt and hr
    if(T){
      dat_plot$pairs_CI <- paste0(dat_plot$HR, " (", dat_plot$HRCILL, " ~ ",dat_plot$HRCIUL, ")")
      dat_plot$p = as.character(dat_plot$p)
      dat_plot$p[c(1,3,4,5,7)] = '<0.001*'
      txt = dat_plot[,c('Varnames','p','pairs_CI')]
      txt <- rbind(c(NA,"pvalue","Hazard ratio"),txt) 
      txt[2,3] = "2.290 (1.500 ~ 3.496)"
      txt[5,3] = "1.462 (1.264 ~ 1.690)"
      txt[7,3] = "1.910 (0.842 ~ 4.332)"
      txt[8,3] = "3.540 (2.363 ~ 5.303)"
      hr <- dat_plot[,c(3:5)] 
      hr <- rbind(rep(NA,ncol(hr)),hr) 
    }
    
    #2.plot
    if(T){
      pdf(file.path(Figure_dir,'step50_forest_uniCOX_TCGA_8.6cm.pdf'),
          width = 3.4, 
          height = 2.8)
      forestplot(txt, hr,
                 graph.pos = 4, 
                 zero = NA, 
                 grid = structure(1.0, gp = gpar(col = "black",lty=2,lwd=2)), 
                 is.summary=c(TRUE,rep(FALSE,7)), 
                 lty.ci = 7, 
                 lwd.ci = 2, 
                 ci.vertices = TRUE, 
                 ci.vertices.height = 0.1, 
                 clip=c(0,6), 
                 fn.ci_norm = fpDrawCircleCI, 
                 boxsize = 0.2, 
                 xlog=FALSE, 
                 xlab = "Hazard ratio", 
                 xticks = c(0,1,2,3,4,5,6), 
                 title = 'Univariate Cox',
                 lwd.xaxis = 1.5, 
                 graphwidth = unit(0.4, "npc"), 
                 colgap = unit(0.04,"npc"), 
                 col=fpColors(box="green", 
                              line="darkblue", 
                              zero='black'), 
                 txt_gp = fpTxtGp(ticks = gpar(cex = 0.7), 
                                  xlab = gpar(cex = 0.7), 
                                  title = gpar(cex = 0.9,fontface = "plain"),
                                  cex = 0.5), 
      )
      dev.off()
    }
  }
}

#4.Multivariate Cox regression  
if(T){
  model <- coxph(Surv(time,status) ~ age+gender+stage+T+N+M+riskScore,
                 data = dat)
  summary(model)
  beta <- coef(model)
  se <- sqrt(diag(vcov(model)))
  HR <- exp(coef(model))
  HRse <- HR * se
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  dat_plot = as.data.frame(tmp[,c('p','HR','HRCILL','HRCIUL')])
  dat_plot = cbind(rownames(tmp),dat_plot)
  colnames(dat_plot)[1] = 'Varnames'
  library(Hmisc)
  dat_plot[,1] = capitalize(dat_plot[,1])
  
  #forestplot
  if(T){
    #1.cpnstruct data txt and hr
    if(T){
      dat_plot$pairs_CI <- paste0(dat_plot$HR, " (", dat_plot$HRCILL, " ~ ",dat_plot$HRCIUL, ")")
      dat_plot$p = as.character(dat_plot$p)
      dat_plot$p[c(4,5,7)] = '<0.001*'
      dat_plot$p[1] = '0.008*'
      txt = dat_plot[,c('Varnames','p','pairs_CI')]
      txt <- rbind(c(NA,"pvalue","Hazard ratio"),txt) 
      txt[5,3] = "1.410 (1.196 ~ 1.662)"
      txt[6,3] = "1.660 (1.310 ~ 2.103)"
      hr <- dat_plot[,c(3:5)] 
      hr <- rbind(rep(NA,ncol(hr)),hr) 
    }
    
    #2.plot
    if(T){
      pdf(file.path(Figure_dir,'step50_forest_multiCOX_TCGA_8.6cm.pdf'),
          width = 3.4, 
          height = 2.8)
      forestplot(txt, hr,
                 graph.pos = 4, 
                 zero = NA, 
                 grid = structure(1.0, gp = gpar(col = "black",lty=2,lwd=2)),
                 is.summary=c(TRUE,rep(FALSE,7)), 
                 lty.ci = 7, 
                 lwd.ci = 2, 
                 ci.vertices = TRUE, 
                 ci.vertices.height = 0.1, 
                 clip=c(0,6), 
                 fn.ci_norm = fpDrawCircleCI, 
                 boxsize = 0.2, 
                 xlog=FALSE, 
                 xlab = "Hazard ratio", 
                 xticks = c(0,1,2,3,4,5,6),  
                 title = 'Multivariate Cox',
                 lwd.xaxis = 1.5, 
                 graphwidth = unit(0.4, "npc"), 
                 colgap = unit(0.04,"npc"), 
                 col=fpColors(box="green", 
                              line="darkblue", 
                              zero='black'), 
                 txt_gp = fpTxtGp(ticks = gpar(cex = 0.7), 
                                  xlab = gpar(cex = 0.7), 
                                  title = gpar(cex = 0.9,fontface = "plain"),
                                  cex = 0.5), 
      )
      dev.off()
    }
  }
}

#5.Construct nomogram
if(T){
  library(survival)
  library(rms)
  library(survivalROC)
  #load the matrix for nomogram construction---dat_nom
  l = load(file = file.path(Rdata_dir,'dat_for_nomogram_construction.Rdata'))
  l
  dd <- datadist(dat_nom)
  options(datadist='dd')
  dd
  
  #Construct surv with function cph()
  coxm <- cph(Surv(time,status) ~ Age+Stage_T+Stage_N+RiskScore,
              data = dat_nom, x=T, y=T, surv = T)
  surv <- Survival(coxm)
  surv1 <- function(x)surv(1,lp=x)
  surv3 <- function(x)surv(3,lp=x)
  surv5 <- function(x)surv(5,lp=x)
  surv10 <- function(x)surv(10,lp=x)
  
  nom<- nomogram(coxm, 
                 fun = list(surv1,surv3,
                            surv5,surv10),
                 fun.at = c(0.1,seq(0.1,0.9,by=0.1),0.9),
                 maxscale = 100,#default maximum point score is 100
                 funlabel = c('1-Year Survival','3-Year Survival',
                              '5-Year Survival','10-Year Survival'),
                 lp=F 
                 )
  Fig_name = file.path(Figure_dir,'step50_nomogram_14_8cm.pdf')
  pdf(Fig_name,
      width = 5.5,
      height = 3.14)
  par(mai = c(0.25,0,0.25,0))
  plot(nom,
       ia.space = 0.2,
       cex.axis = 0.4,#character size for tick mark labels
       cex.var = 0.6, #	character size for axis titles (variable names)
       xfrac = 0.18,#fraction of horizontal plot to set aside for axis titles
       tcl=-0.15,#length of tick marks in nomogram
       lmgp = 0.05 #spacing between numeric axis labels and axis
  )
  dev.off()

}

#6.Internal validation of the nomogram（2 indicators）：
# 1）discrimination：C-index
# 2）calibration：calibration curve
if(T){
  #1）Calculate C-index
  rcorrcens(Surv(time, status) ~ predict(coxm), data = dat_nom)

  #Calculate AUC
  # 1-year AUC
  if(T){
    f1 <- cph(Surv(time, status) ~ Age+Stage_T+Stage_N+RiskScore,
              data = dat_nom,
              x=T, y=T, surv=T, time.inc=1)
    cutoff = 1
    roc_1 = survivalROC(Stime=dat_nom$time, 
                        status=dat_nom$status, #1==dead，0==alive
                        marker = predict(f1),
                        predict.time =cutoff, #predict.time：time point of the ROC curve
                        method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_1)
  }
  
  # 3-year AUC
  if(T){
    f3 <- cph(Surv(time, status) ~ Age+Stage_T+Stage_N+RiskScore,
              data = dat_nom,
              x=T, y=T, surv=T, time.inc=3)
    cutoff = 3
    roc_3 = survivalROC(Stime=dat_nom$time, 
                        status=dat_nom$status, 
                        marker = predict(f3),
                        predict.time =cutoff, #predict.time：time point of the ROC curve
                        method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_3)
  }
  
  # 5-year AUC
  if(T){
    f5 <- cph(Surv(time, status) ~ Age+Stage_T+Stage_N+RiskScore,
              data = dat_nom,
              x=T, y=T, surv=T, time.inc=5)
    cutoff = 5
    roc_5 = survivalROC(Stime=dat_nom$time, 
                        status=dat_nom$status, 
                        marker = predict(f5),
                        predict.time =cutoff, #predict.time：time point of the ROC curve
                        method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_5)
  }
  # 10-year AUC
  if(T){
    f10 <- cph(Surv(time, status) ~ Age+Stage_T+Stage_N+RiskScore,
              data = dat_nom,
              x=T, y=T, surv=T, time.inc=10)
    cutoff = 10
    roc_10 = survivalROC(Stime=dat_nom$time, 
                        status=dat_nom$status, 
                        marker = predict(f10),
                        predict.time =cutoff, #predict.time：time point of the ROC curve
                        method="KM" #Method for fitting joint distribution of (marker,t)
    )
    str(roc_10)
  }
  #plot
  if(T){
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
    
    roc = rbind(roc1,roc3,roc5,roc10)
    if(T){
      p <- ggplot(roc, aes(x=FP, y=TP, colour=group,group=group)) + 
        geom_line(size=0.5)+
        geom_abline(intercept = 0,slope = 1,lty=2,lwd=0.4,colour="black")+
        xlim(0,1.03) + ylim(0,1.03) +
        labs(title="Training cohort",
             x="False positive rate",
             y="True positive rate") +
        scale_colour_manual(values = c('red','blue','green','orange'),
                            limits = c('1','3','5','10'),
                            labels = c(paste("1 year AUC: ",round(roc_1$AUC,3)),
                                       paste("3 year AUC: ",round(roc_3$AUC,3)),
                                       paste("5 year AUC: ",round(roc_5$AUC,3)),
                                       paste("10 year AUC: ",round(roc_10$AUC,3)))
        ) 
      p_roc <- p + theme(
        panel.background = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.grid=element_blank(),
        plot.title    = element_text(color = 'black', size   = 8.5, face = 'bold',hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
        axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
        legend.position=c(0.73,0.18),
        legend.key.size = unit(8, "pt"),
        legend.title  = element_blank(),
        legend.text   = element_text(color = 'black', size   = 6),
        axis.line = element_blank(), 
        panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
      )
      print(p_roc)
      save(p_roc, file = file.path(Rdata_dir,'step50_nomogram_roc.Rdata'))
      ggsave(filename = file.path(Figure_dir,'step50_nomogram_ROC_1-3-5-10year_6-6cm.pdf'),
             plot = last_plot(), width = 6, height = 6, units = 'cm', dpi = 600)
    }
  }
  
  #2）calibration curve
  if(T){
    if(T){
      f1 <- cph(Surv(time, status) ~ Age+Stage_T+Stage_N+RiskScore,
                data = dat_nom,
                x=T, y=T, surv=T, time.inc=1)
      f3 <- cph(Surv(time, status) ~ Age+Stage_T+Stage_N+RiskScore,
                data = dat_nom,
                x=T, y=T, surv=T, time.inc=3)
      f5 <- cph(Surv(time, status) ~ Age+Stage_T+Stage_N+RiskScore,
                data = dat_nom,
                x=T, y=T, surv=T, time.inc=5)
      f10 <- cph(Surv(time, status) ~ Age+Stage_T+Stage_N+RiskScore,
                 data = dat_nom,
                 x=T, y=T, surv=T, time.inc=10)
      cal1 <- calibrate(f1, cmethod="KM", method="boot", 
                        u=1, # u represents the predicted time point, which should be consistent with time.inc defined in above model  
                        m=112, # m should be determined according to sample size. 
                        #The standard curve generally divides all samples into 3 groups (3 points shown in the figure).  
                        #m represents the sample size of each group, so m*3 should be equal to or approximately equal to the sample size
                        B=336 
      )
      cal3 <- calibrate(f3, cmethod="KM", method="boot", 
                        u=3, 
                        m=112, 
                        B=336 
      )
      cal5 <- calibrate(f5, cmethod="KM", method="boot", 
                        u=5, 
                        m=112, 
                        B=336 
      )
      cal10 <- calibrate(f10, cmethod="KM", method="boot", 
                         u=10, 
                         m=112,
                         B=500 
      )
    }
    
    # Plot
    if(T){
      Fig_name = file.path(Figure_dir,'step09_nom_calibration_6_6cm.pdf')
      pdf(Fig_name,
          width = 4,
          height = 4)
      #mai = c(bottom, left, top, right) in inches.
      par(mai = c(0.8,0.8,0.05,0.1))
      plot(cal1,lwd=1.2,lty=1,
           cex.axis=0.7, 
           cex.lab=0.8, 
           errbar.col='red',
           xlim=c(0,1.05),
           ylim=c(0,1.05),
           xlab="Predicted Survival Probability",
           ylab="Actual Survival Probability",
           col='red',
           subtitles=F 
      )
      plot(cal3,lwd=1,lty=1,
           errbar.col='blue',col='blue',
           add=T #set to TRUE to add the calibration plot to an existing plot
      )
      plot(cal5,lwd=1,lty=1,errbar.col='green',col='green',add=T)
      plot(cal10,lwd=1,lty=1,errbar.col='orange',col='orange',add=T)
      text(0.86,0.03,'10-year Survival',col='orange',cex =0.7)
      text(0.86,0.1,'5-year Survival',col='green',cex =0.7)
      text(0.86,0.17,'3-year Survival',col='blue',cex =0.7)
      text(0.86,0.24,'1-year Survival',col='red',cex =0.7)
      dev.off()
    }
  }
}

#7.External validation of the nomogram with GEO dataset
if(T){
  GSEName = 'GSE98394'
  # step01.load clinical data and riskScore of GEO dataset
  if(T){
    load(file = file.path(Rdata_dir,paste(GSEName,'RData',sep=".")))
    a=as.data.frame(t(exprSet))
    ggplot(a,aes(x=a[,2]))+geom_density()
    
    identical(rownames(meta),colnames(exprSet))
    phe = na.omit(meta[,c(11,12,4,7,8)])
    colnames(phe) = c('time','status','age','T','N')
    phe$time = as.numeric(phe$time)
    phe$time = phe$time/12 
    phe = na.omit(phe) 
    dim(phe) 
    exprSet = exprSet[,rownames(phe)]
    dim(exprSet) 
    #Calculate RS for patients
    if(T){
      load(file = file.path(Rdata_dir,'step50_LASSO_modeleGenes_coef.RData'))
      LASSOgenes = rownames(coef_out3)
      exprSet_t = t(exprSet)
      LASSOgenes[3] = 'RASSF3'
      expr_LASSOgenes = exprSet_t[,LASSOgenes]
      phe = cbind(phe,expr_LASSOgenes)
      
      phe$riskScore = coef_out3[1,1]*phe[,6] + coef_out3[2,1]*phe[,7] +
        coef_out3[3,1]*phe[,8] + coef_out3[4,1]*phe[,9] +
        coef_out3[5,1]*phe[,10] + coef_out3[6,1]*phe[,11] + coef_out3[7,1]*phe[,12]
      phe = phe[order(phe$riskScore,decreasing = FALSE),]
    }
    
    #Collect data for medel validation
    dat_geo = phe
    dat_geo$Age = ifelse(dat_geo$age<75,'<75','>=75')
    dat_geo$Age = factor(dat_geo$Age,ordered=F,levels = c('<75','>=75'))
    dat_geo$Stage_T = factor(substr(dat_geo$T,1,2),ordered=F,levels = c('T0','T1','T2','T3','T4'))
    stage_N = substr(dat_geo$N,1,2)
    dat_geo$Stage_N = ifelse(stage_N=='Nx','N1',stage_N)  
    colnames(dat_geo)[13] = 'RiskScore'
    dat_geo = dat_geo[,c('time','status','Age','Stage_T','Stage_N','RiskScore')]
  }
  
  # step02.C-index
  if(T){
    fev <- cph(Surv(time, status) ~ predict(coxm, newdata=dat_geo), 
               x=T, y=T, surv=T, data=dat_geo)  
    rcorrcens(Surv(time, status) ~ predict(coxm, newdata=dat_geo), data = dat_geo)
    #Another way to calculate C-index
    if(T){
      v = validate(fev, method="boot", B=1000, dxy=T)
      Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
      orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"]
      
      # The c-statistic according to Dxy=2(c-0.5)
      bias_corrected_c_index  <- abs(Dxy)/2+0.5
      orig_c_index <- abs(orig_Dxy)/2+0.5
      
      bias_corrected_c_index
      orig_c_index
    }
    if(T){
      library(survcomp)
      cindex <- concordance.index(predict(fev),
                                  surv.time = dat_geo$time, 
                                  surv.event = dat_geo$status,
                                  method = "noether")
      cindex$c.index; cindex$lower; cindex$upper
    }
    
    #Calculate AUC
    # 1-year AUC
    if(T){
      fev1 <- cph(Surv(time, status) ~ predict(coxm, newdata=dat_geo),
                  data = dat_geo,
                  x=T, y=T, surv=T, time.inc=1)
      cutoff = 1
      roc_1 = survivalROC(Stime=dat_geo$time, 
                          status=dat_geo$status, #1==dead，0==alive
                          marker = predict(fev1),
                          predict.time =cutoff, #predict.time：time point of the ROC curve
                          method="KM" #Method for fitting joint distribution of (marker,t)
      )
      str(roc_1)
    }
    # 3-year AUC
    if(T){
      fev3 <- cph(Surv(time, status) ~ predict(coxm, newdata=dat_geo),
                  data = dat_geo,
                  x=T, y=T, surv=T, time.inc=3)
      cutoff = 3
      roc_3 = survivalROC(Stime=dat_geo$time, 
                          status=dat_geo$status, 
                          marker = predict(fev3),
                          predict.time =cutoff, #predict.time：time point of the ROC curve
                          method="KM" #Method for fitting joint distribution of (marker,t)
      )
      str(roc_3)
    }
    # 5-year AUC
    if(T){
      fev5 <- cph(Surv(time, status) ~ predict(coxm, newdata=dat_geo),
                  data = dat_geo,
                  x=T, y=T, surv=T, time.inc=5)
      cutoff = 5
      roc_5 = survivalROC(Stime=dat_geo$time, 
                          status=dat_geo$status, 
                          marker = predict(fev5),
                          predict.time =cutoff, #predict.time：time point of the ROC curve
                          method="KM" #Method for fitting joint distribution of (marker,t)
      )
      str(roc_5)
    }
    # 10-year AUC
    if(T){
      fev10 <- cph(Surv(time, status) ~ predict(coxm, newdata=dat_geo),
                  data = dat_geo,
                  x=T, y=T, surv=T, time.inc=10)
      cutoff = 10
      roc_10 = survivalROC(Stime=dat_geo$time, 
                          status=dat_geo$status, 
                          marker = predict(fev10),
                          predict.time =cutoff, #predict.time：time point of the ROC curve
                          method="KM" #Method for fitting joint distribution of (marker,t)
      )
      str(roc_10)
    }
    #Plot
    if(T){
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
      
      roc = rbind(roc3,roc5,roc10)
      if(T){
        p <- ggplot(roc, aes(x=FP, y=TP, colour=group,group=group)) + 
          geom_line(size=0.5)+
          geom_abline(intercept = 0,slope = 1,lty=2,lwd=0.4,colour="black")+
          xlim(0,1.03) + ylim(0,1.03) +
          labs(title="Validation cohort",
               x="False positive rate",
               y="True positive rate") +
          scale_colour_manual(values = c('blue','green','orange'),
                              limits = c('3','5','10'),
                              labels = c(paste("3 year AUC: ",round(roc_3$AUC,3)),
                                         paste("5 year AUC: ",round(roc_5$AUC,3)),
                                         paste("10 year AUC: ",round(roc_10$AUC,3)))
          ) 
        p_roc <- p + theme(
          panel.background = element_rect(fill = 'transparent'), 
          legend.background = element_rect(fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent'),
          panel.grid=element_blank(),
          plot.title    = element_text(color = 'black', size   = 8.5, face='bold',hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
          plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
          axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
          axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
          axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
          axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
          legend.position=c(0.73,0.15),
          legend.key.size = unit(8, "pt"),
          legend.title  = element_blank(),
          legend.text   = element_text(color = 'black', size   = 6),
          axis.line = element_blank(), 
          panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
        )
        print(p_roc)
        save(p_roc, file = file.path(Rdata_dir,'step50_nomogram_roc_validation-GEO.Rdata'))
        ggsave(filename = file.path(Figure_dir,'step50_nomogram_ROC_validation-GEO_3-5-10year_6-6cm.pdf'),
               plot = last_plot(), width = 6, height = 6, units = 'cm', dpi = 600)
      }
    }
  }
  
  # step03.calibrate curve 
  if(T){
    if(T){
      #Only 1 Sample with OS<1, so fev1 can't be calculated
      fev3 <- cph(Surv(time, status) ~ predict(coxm, newdata=dat_geo),
                data = dat_geo,
                x=T, y=T, surv=T, time.inc=3)
      fev5 <- cph(Surv(time, status) ~ predict(coxm, newdata=dat_geo),
                  data = dat_geo,
                  x=T, y=T, surv=T, time.inc=5)
      fev10 <- cph(Surv(time, status) ~ predict(coxm, newdata=dat_geo),
                  data = dat_geo,
                  x=T, y=T, surv=T, time.inc=10)
      calev3 <- calibrate(fev3, cmethod="KM", method="boot", 
                          u=3, m=14, B=1000)
      calev5 <- calibrate(fev5, cmethod="KM", method="boot", 
                          u=5, m=14, B=1000)
      calev10 <- calibrate(fev10, cmethod="KM", method="boot", 
                          u=10, m=14, B=1000)
    }
    ## plot
    if(T){
      Fig_name = file.path(Figure_dir,'step50_nom_calibration_validation-GEO.pdf')
      pdf(Fig_name,
          width = 4,
          height = 4)
      #mai = c(bottom, left, top, right) in inches.
      par(mai = c(0.8,0.8,0.05,0.1))
      plot(calev3,lwd=1.2,lty=1,
           cex.axis=0.7, 
           cex.lab=0.8, 
           errbar.col='red',
           xlim=c(0,1.05),
           ylim=c(0,1.05),
           xlab="Predicted Survival Probability",
           ylab="Actual Survival Probability",
           col='blue',
           subtitles=F 
      )
      plot(calev5,lwd=1.2,lty=1,errbar.col='green',col='green',add=T)
      plot(calev10,lwd=1.2,lty=1,errbar.col='orange',col='orange',add=T)
      text(0.86,0.03,'10-year Survival',col='orange',cex =0.7)
      text(0.86,0.1,'5-year Survival',col='green',cex =0.7)
      text(0.86,0.17,'3-year Survival',col='blue',cex =0.7)
      dev.off()
    }
}
}




