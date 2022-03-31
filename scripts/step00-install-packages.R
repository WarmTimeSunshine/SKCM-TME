### install packages

rm(list = ls()) 
options(warn =-1) 
options()$repos  
options()$BioC_mirror 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
# Specified download path with BioCManager
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
options()$repos 
options()$BioC_mirror


# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 

if(!require("patchwork")) BiocManager::install("patchwork",ask = F,update = F)
if(!require("cowplot")) BiocManager::install("cowplot",ask = F,update = F)
if(!require("KEGG.db")) BiocManager::install("KEGG.db",ask = F,update = F)
if(!require("GSEABase")) BiocManager::install("GSEABase",ask = F,update = F)
if(!require("GSVA")) BiocManager::install("GSVA",ask = F,update = F)
if(!require("clusterProfiler")) BiocManager::install("clusterProfiler",ask = F,update = F)
if(!require("GEOquery")) BiocManager::install("GEOquery",ask = F,update = F)
if(!require("limma")) BiocManager::install("limma",ask = F,update = F)
if(!require("DESeq2")) BiocManager::install("DESeq2",ask = F,update = F)
if(!require("edgeR")) BiocManager::install("edgeR",ask = F,update = F)
if(!require("airway")) BiocManager::install("airway",ask = F,update = F)
if(!require("impute")) BiocManager::install("impute",ask = F,update = F)
if(!require("genefu")) BiocManager::install("genefu",ask = F,update = F)
if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db",ask = F,update = F)
if(!require("hgu133plus2.db")) BiocManager::install("hgu133plus2.db",ask = F,update = F)
if(!require("ConsensusClusterPlus")) BiocManager::install("ConsensusClusterPlus",ask = F,update = F)
if(!require("biomaRt"))BiocManager::install("biomaRt",ask = F,update = F)
if(!require("annotate"))BiocManager::install("annotate",ask = F,update = F)
if(! require("maftools")) BiocManager::install("maftools",ask = F,update = F)
if(! require("genefilter")) BiocManager::install("genefilter",ask = F,update = F)
if(! require("CLL")) BiocManager::install("CLL",ask = F,update = F)
if(! require("RTCGA")) BiocManager::install("RTCGA",ask = F,update = F)
if(! require("RTCGA.clinical")) BiocManager::install("RTCGA.clinical",ask = F,update = F)
# https://bioconductor.org/packages/3.6/data/experiment/src/contrib/RTCGA.clinical_20151101.8.0.tar.gzn
if(! require("RTCGA.miRNASeq")) BiocManager::install("RTCGA.miRNASeq",ask = F,update = F)
if(! require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks",ask = F,update = F)

if(!require("WGCNA")) install.packages("WGCNA",update = F,ask = F)
if(!require("FactoMineR")) install.packages("FactoMineR",update = F,ask = F)
if(!require("factoextra")) install.packages("factoextra",update = F,ask = F)
if(!require("ggplot2")) install.packages("ggplot2",update = F,ask = F)
if(!require("pheatmap")) install.packages("pheatmap",update = F,ask = F)
if(!require("ggpubr")) install.packages("ggpubr",update = F,ask = F)
if(!require("glmnet")) install.packages("glmnet",update = F,ask = F)
if(!require("VennDiagram")) install.packages("VennDiagram",update = F,ask = F)
if(!require("randomForest")) install.packages("randomForest",update = F,ask = F)
if(!require("forestplot")) install.packages("forestplot",update = F,ask = F)
install.packages(c("devtools","reshape2","ggfortify","stringr",
                   "survival","survminer","lars","timeROC",
                   "ROCR","Hmisc","caret","ggstatsplot","tableone"))

library("FactoMineR")
library("factoextra")
library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(genefu)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)
library(devtools) 

