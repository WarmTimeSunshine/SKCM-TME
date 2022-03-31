#LASSO regression analysis to screen genes for model construction

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'

#STEP 01.load expression data, clinical information,and gene names for analysis
if(T){
  workflow_type = "HTSeq - FPKM"
  exprSym_f=file.path(Rdata_dir,
                      paste(paste("TCGA",request_cancer,workflow_type,'exprSym_meta_groupT-N.RData',sep="-")))
  load(file = exprSym_f)
  exprSym_FPKM = na.omit(exprSym_FPKM)
  dim(exprSym_FPKM)
  dim(meta)
  table(group_list_TN)
  
  GSEName = 'GSE65904'
  surv_f=file.path(Rdata_dir,
                   paste('survival_results_cox_TCGA.RData',GSEName,sep="-"))
  load(file = surv_f)
  sigGenes_cox = sigGenes_cox_GEO
  
  #Procesing of exprSym_FPKM
  exprDat = exprSym_FPKM[rownames(sigGenes_cox),group_list_TN=="Tumor"]
  exprDat = exprDat[,!duplicated(substr(colnames(exprDat),1,12))]
  dim(exprDat)
  exprDat = as.data.frame(exprDat)
  
  #Procesing of meta
  sampleID = substr(colnames(exprDat),1,12)
  table(table(sampleID)) 
  clin_info = meta[meta$submitter_id %in% sampleID,
                   c('submitter_id','vital_status','days_to_death','days_to_last_follow_up')]
  clin_info = as.data.frame(clin_info) 
  clin_info$time = ifelse(clin_info$vital_status=="Alive",
                          clin_info$days_to_last_follow_up, clin_info$days_to_death)
  clin_info$status=ifelse(clin_info$vital_status=="Dead",1,0)
  rownames(clin_info) = clin_info$submitter_id
  phe = na.omit(clin_info[,c('time','status')])
  phe = phe[phe$time>30,]
  phe$time = phe$time/365 
  #phe includes 2 columns: time and event
  
  exprSet = exprDat[,match(rownames(phe),substr(colnames(exprDat),1,12))]
  exprSet = log2(exprSet+1)
  dim(exprSet)
  save(exprSet,phe, file = file.path(Rdata_dir,'step50_LASSO_input_data.RData'))
}

#STEP 02.lasso regression analysis
if(T){
  library(lars) 
  library(glmnet) 
  library(ggplot2)
  
  identical(substr(colnames(exprSet),1,12),rownames(phe))
  x=as.matrix(t(exprSet))
  rownames(x)<-NULL
  phe$time = as.double(phe$time)
  phe$status = as.double(phe$status)
  y=as.matrix(survival::Surv(phe$time,phe$status))
  
  #Cross-validation is used to find the best modeling location
  if(T){
    set.seed(12345)
    cv_fit <- cv.glmnet(x, y, family='cox', nfolds = 10, alpha = 1)
    pdf(file.path(Figure_dir,'step50_LASSO_cross_validation.pdf'),
        width = 6,  
        height = 6  
    )
    par(mai=c(1,1,0.5,0.25))#c(bottom, left, top, right)
    plot(cv_fit,
         cex=0.01, 
         xlab = 'Log(Lambda)',
         cex.axis=1.4, 
         cex.lab=1.6 
    )
    dev.off()
    # The two dashed lines indicate two special values of λ,cv_fit$lambda.1se is the Optimal modeling location
    c(cv_fit$lambda.min,cv_fit$lambda.1se) 

    #coef is the coefficient of each gene in the model
    coef_lambda_min = coef(cv_fit, s=cv_fit$lambda.min)
    coef_lambda_min
    coef_out = coef_lambda_min[which(coef_lambda_min!=0),]
    coef_out2 = matrix(coef_out,length(coef_out),1)
    rownames(coef_out2)=names(coef_out)
    colnames(coef_out2)=c('coef')
    coef_out2
    coef_out3 = as.matrix(coef_out2[order(abs(coef_out2),decreasing = TRUE),])
    colnames(coef_out3)=c('coef')  
    save(coef_out3, file = file.path(Rdata_dir,'step50_LASSO_modeleGenes_coef.RData'))
    
    model_lasso <- glmnet(x, y, family="cox", alpha = 1)
    print(model_lasso)
    # The column %Dev Represents the proportion of residuals explained by the model, 
    # equal to R^2 for linear models  
    # The higher the value, the better the model performed
    pdf(file.path(Figure_dir,'step50_LASSO_coefficient_profiles.pdf'),
        width = 6, 
        height = 6 
    )
    par(mai=c(1,1,0.5,0.25))#c(bottom, left, top, right)
    plot(model_lasso, xvar="lambda", label=F,
         xlab = 'Log(Lambda)',
         cex.axis=1.4, 
         cex.lab=1.6 
    )
    abline(v=log(cv_fit$lambda.min, base = exp(1)),
           lwd=1.8, lty = 3,
           col="black")
    dev.off()
  }
}

#STEP 03.Group samples based on risk score
# Plot RS distribution，heatmap，KMcurve
if(T){
  load(file = file.path(Rdata_dir,'step50_LASSO_modeleGenes_coef.RData'))
  load(file = file.path(Rdata_dir,'step50_LASSO_input_data.RData'))
  dim(exprSet)
  dim(phe)
  identical(rownames(phe),substr(colnames(exprSet),1,12))
  genes = rownames(coef_out3)
  exprSet_t = t(exprSet)
  expr_lassoGenes = exprSet_t[,genes]
  phe = cbind(phe,expr_lassoGenes)
  phe$riskScore = coef_out3[1,1]*phe[,3] + coef_out3[2,1]*phe[,4] +
    coef_out3[3,1]*phe[,5] + coef_out3[4,1]*phe[,6] +
    coef_out3[5,1]*phe[,7] + coef_out3[6,1]*phe[,8] + coef_out3[7,1]*phe[,9] 
  phe$RS_group = ifelse(phe$riskScore>median(phe$riskScore),'high','low')
  phe = phe[order(phe$riskScore,decreasing = FALSE),]
  phe$samNum = 1:nrow(phe)
  phe$event = ifelse(phe$status==0,'Alive','Dead')
  save(phe, file = file.path(Rdata_dir,'step50_LASSO_KM_input.RData'))
  
  a=phe
  
  #1)RS distribution
  if(T){
    library(ggplot2)
    plot.point=ggplot(a,aes(x=samNum,y=riskScore,color=RS_group)) +
      geom_point(size=1) +
      scale_colour_manual(values = c('green','red'),limits = c('low','high'),
                          labels = c('Low risk','High risk'))+ 
      xlab("Patients (increasing risk score)") +
      ylab("Risk score") +
      geom_hline(yintercept = median(a$riskScore),lty=2,lwd=0.4,alpha=0.8) +
      geom_vline(xintercept = nrow(a[a$RS_group=='low',]),lty=2,lwd=0.4,alpha=0.8)+
      theme(
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
        legend.position=c(0.13,0.88),
        legend.key.size = unit(7.5, "pt"),
        legend.title  = element_blank(),
        legend.text   = element_text(color = 'black', size   = 6),
        axis.line = element_blank(), 
        panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
      )
    print(plot.point)
    ggsave(filename = file.path(Figure_dir,'step50_LASSO_model_riskScore_distribution_8-4cm.pdf'),
           plot = last_plot(), width = 8, height = 4, units = 'cm', dpi = 600)
  }
  
  #2)Dot plot of survival status
  if(T){
    plot.sur=ggplot(phe,aes(x=samNum,y=time))+
      geom_point(aes(col=event),size=0.65) +
      xlab("Patients (increasing risk score)") +
      ylab("Survival time (years)") +
      theme(
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
        legend.position=c(0.9,0.92),
        legend.key.size = unit(7.5, "pt"),
        legend.title  = element_blank(),
        legend.text   = element_text(color = 'black', size   = 6),
        axis.line = element_blank(), 
        panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
      )
    print(plot.sur)
    ggsave(filename = file.path(Figure_dir,'step50_LASSO_model_survival_time_8-4cm.pdf'),
           plot = last_plot(), width = 8, height = 4, units = 'cm', dpi = 600)
  }
  
  #3)heatmap 
  if(T){
    mycolors <- colorRampPalette(c("green", "black", "red"), bias = 1.2)(100)
    exp_dat = phe[,genes]
    tmp=t(scale(exp_dat,center = TRUE, scale = TRUE))
    tmp[tmp > 1] = 1
    tmp[tmp < -1] = -1
    
    library(pheatmap)
    plot.h=pheatmap(tmp,col= mycolors,
                    show_colnames = F,
                    show_rownames = F,
                    cluster_cols = F,
                    cluster_rows = F,
                    fontsize = 6)
    ggsave(filename = file.path(Figure_dir,'step50_LASSO_model_heatmap_10-8cm.pdf'),
           plot = plot.h, width = 10, height = 8, units = 'cm', dpi = 600)
  }
}






