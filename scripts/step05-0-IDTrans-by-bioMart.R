#Introduction:This script converts rownames in the gene expression matrix from ensemble ID to Symbol ID. 
#The matrix obtained after transformation is named exprSym_FPKM.
#exprSym_FPKM,clinical data (meta) and grouping data (tumor,normal)   
#are save together named 'TCGA-COAD-FPKM_exprSym_meta_group.Rdata'  

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'

# load mRNA and clinical data downloaded from GDC
tmpfile = file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,data_type,workflow_type,'getTCGAData',sep="-"),
                          'RData',sep="."))
load(file = tmpfile)
dim(expr)
dim(meta)

#ID transform with bioMart----
library(biomaRt)
options(RCurlOptions = list(proxy="uscache.kcc.com:80",proxyuserpwd="------:-------"))
martFile = file.path(Rdata_dir,"bioMart-hsapiens_ensembl.Rdata")
if(!file.exists(martFile)){
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  save(mart, file = martFile)
}
load(file = martFile)
gene_symbols<- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                     filters= "ensembl_gene_id", 
                     values = rownames(expr), 
                     mart = mart)
table(table(gene_symbols$hgnc_symbol))
table(table(gene_symbols$ensembl_gene_id))

#After trasformation----
r <- with(gene_symbols, which(hgnc_symbol=="", arr.ind=TRUE))
gene_symbols <- gene_symbols[-r, ]

duplicate_ensembl <- gene_symbols$ensembl_gene_id[duplicated(gene_symbols$ensembl_gene_id)]
gene_symbols <- gene_symbols[!(gene_symbols$ensembl_gene_id %in% duplicate_ensembl), ]
duplicate_symbol <- gene_symbols$hgnc_symbol[duplicated(gene_symbols$hgnc_symbol)]
gene_symbols <- gene_symbols[!(gene_symbols$hgnc_symbol %in% duplicate_symbol), ]

expr_toMerge <- cbind(rownames(expr), expr)
exprSym_FPKM = merge(gene_symbols, expr_toMerge, by.x="ensembl_gene_id", by.y='V1')
rownames(exprSym_FPKM) <- exprSym_FPKM$hgnc_symbol
exprSym_FPKM <- as.matrix(exprSym_FPKM[,-c(1,2)])
exprSym_FPKM = apply(exprSym_FPKM, c(1,2), as.numeric)
dim(exprSym_FPKM) 
save(exprSym_FPKM, file = file.path(Rdata_dir,"TCGA-SKCM_expr_FPKM.Rdata"))

#Load counts data and meta
exprSym_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,'exprSym_meta_groupT-N.RData',sep="-")))
load(file = exprSym_f)
table(group_list_TN)
dim(exprSym)

zero_num <- apply(exprSym_FPKM, 1, function(gene){
  n = length(which((gene==0)))
  return(n)
})
length(zero_num)
table(zero_num)
expr_FPKM = exprSym_FPKM[zero_num<(ncol(exprSym_FPKM)*0.2),]
expr_FPKM = na.omit(expr_FPKM)
dim(expr_FPKM)
identical(rownames(exprSym),rownames(expr_FPKM))

expr_FPKM = expr_FPKM[,(substr(colnames(expr_FPKM),1,12) %in% meta$submitter_id)]
dim(expr_FPKM)
exprSym_FPKM = expr_FPKM

#To get group data（tumor,normal）----
table(substring(colnames(exprSym_FPKM),14,15)) 
sample_num <- as.numeric(substring(colnames(exprSym_FPKM),14,15)) 
table(sample_num)
group_list_TN=ifelse(sample_num %in% 1:9,"Tumor","Normal")
table(group_list_TN)

#Management of exprSym_FPKM and group_list_TN
if(T){
  exprSym_tumor = exprSym_FPKM[,group_list_TN=='Tumor']
  dim(exprSym_tumor) 
  exprSym_tumor = exprSym_tumor[,!duplicated(substr(colnames(exprSym_tumor),1,15))]
  dim(exprSym_tumor) 
  exprSym_FPKM = exprSym_tumor
  dim(exprSym_FPKM) 
}

#Save data----
exprSym_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,workflow_type,'exprSym_meta_groupT-N.RData',sep="-")))
save(exprSym_FPKM, meta, group_list_TN, file = exprSym_f)


