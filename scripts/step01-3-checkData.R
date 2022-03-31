#Examine data before DEG

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'

#load data:exprSym，meta (clinical information)，group_list_TN
if(T){
  exprSym_f=file.path(Rdata_dir,
                      paste(paste("TCGA",request_cancer,'exprSym_meta_groupT-N.RData',sep="-")))
  load(file = exprSym_f)
  dim(exprSym)
  exprSet=na.omit(exprSym)
  dim(exprSet)
  table(group_list_TN)
}

#STEP 1:Examine common genes（such as GAPDH、ACTB）
if(T){
  boxplot(exprSet[1,])
  boxplot(exprSet[2,])
  boxplot(exprSet[3,])
  boxplot(exprSet['GAPDH',])
  boxplot(exprSet['ACTB',])
}

#STEP 2:distribution of gene expression
if(T){
  library(reshape2)
  exprSet_L=melt(exprSet)
  colnames(exprSet_L)=c('symbol','sample','value')
  exprSet_L$group=rep(group_list_TN,each=nrow(exprSet))
  head(exprSet_L)
  library(ggplot2)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
  print(p)
  p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density() 
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p)
}

#STEP 3:Examine the group information,such as hclust,PCA
if(T){
  # hclust
  exprSet_hclust=exprSet
  colnames(exprSet_hclust)=paste(group_list_TN,1:ncol(exprSet),sep='')
  plot(hclust(dist(t(exprSet_hclust))))
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet_hclust)))
  par(mar=c(5,5,5,10)) 
  plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)
  
  # PCA
  #BiocManager::install('ggfortify')
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list_TN 
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
  library("FactoMineR")
  library("factoextra") 
  df=as.data.frame(t(exprSet))
  dat.pca <- PCA(df, graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list_TN, # color by groups
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  ggsave('all_samples_PCA.png')
}