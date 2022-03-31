#Introduction:This script converts rownames in the gene expression matrix from ensemble ID to Symbol ID. 
#The matrix obtained after transformation is named exprSym.
#exprSym,clinical data (meta) and grouping data (tumor,normal)   
#are save together named 'TCGA-SKCM-exprSym_meta_group.Rdata'  

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'

# load mRNA and clinical data
tmpfile = file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,data_type,workflow_type,'getTCGAData',sep="-"),
                          'RData',sep="."))
load(file = tmpfile)
dim(expr)
dim(meta)

#ID trasformation----
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

#After trasformation
r <- with(gene_symbols, which(hgnc_symbol=="", arr.ind=TRUE))
gene_symbols <- gene_symbols[-r, ]

duplicate_ensembl <- gene_symbols$ensembl_gene_id[duplicated(gene_symbols$ensembl_gene_id)]
gene_symbols <- gene_symbols[!(gene_symbols$ensembl_gene_id %in% duplicate_ensembl), ]
duplicate_symbol <- gene_symbols$hgnc_symbol[duplicated(gene_symbols$hgnc_symbol)]
gene_symbols <- gene_symbols[!(gene_symbols$hgnc_symbol %in% duplicate_symbol), ]

expr_toMerge <- cbind(rownames(expr), expr)
exprSym = merge(gene_symbols, expr_toMerge, by.x="ensembl_gene_id", by.y='V1')
rownames(exprSym) <- exprSym$hgnc_symbol
exprSym <- as.matrix(exprSym[,-c(1,2)])
exprSym = apply(exprSym, c(1,2), as.numeric)
dim(exprSym) 

zero_num <- apply(exprSym, 1, function(gene){
  n = length(which((gene==0)))
  return(n)
})
length(zero_num)
table(zero_num)
exprSym = exprSym[zero_num<(ncol(exprSym)*0.2),]
exprSym = na.omit(exprSym)
dim(exprSym)

if(T){
  clin_info = meta[,c('submitter_id','ajcc_pathologic_stage',
                      'ajcc_pathologic_t','ajcc_pathologic_n','ajcc_pathologic_m',
                      'vital_status','days_to_death','days_to_last_follow_up')]
  #Overall survical:days_to_last_follow_up if alive，days_to_death if dead
  clin_info$time = ifelse(clin_info$vital_status=="Alive",
                          clin_info$days_to_last_follow_up, clin_info$days_to_death)
  #Dead==1，Alive==0
  clin_info$status=ifelse(clin_info$vital_status=="Dead",1,0)
  clin_info = clin_info[,-c(7,8)]
  clin_info = na.omit(clin_info)
  clin_info = clin_info[clin_info$ajcc_pathologic_stage!='Not Reported',]
  clin_info = clin_info[clin_info$ajcc_pathologic_stage!='Stage 0',]
  clin_info = clin_info[clin_info$ajcc_pathologic_t!='TX',]
  clin_info = clin_info[clin_info$ajcc_pathologic_t!='Tis',]
  clin_info = clin_info[clin_info$ajcc_pathologic_n!='NX',]
  clin_info = clin_info[clin_info$time>30,]
  median(clin_info$time)
  
  meta = meta[meta$submitter_id %in% clin_info$submitter_id,]
  exprSym = exprSym[,substr(colnames(exprSym),1,12) %in% clin_info$submitter_id]
  dim(exprSym)
}

#group information（tumor,normal）----
table(substring(colnames(exprSym),14,15)) 
sample_num <- as.numeric(substring(colnames(exprSym),14,15)) 
table(sample_num)
group_list_TN=ifelse(sample_num %in% 1:9,"Tumor","Normal")
table(group_list_TN) 

#Management of exprSym and group_list_TN
if(T){
  exprSym_tumor = exprSym[,group_list_TN=='Tumor']
  dim(exprSym_tumor) 
  exprSym_tumor = exprSym_tumor[,!duplicated(substr(colnames(exprSym_tumor),1,15))]
  dim(exprSym_tumor) 
  dim(exprSym) 
}

#Save data----
exprSym_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,'exprSym_meta_groupT-N.RData',sep="-")))
save(exprSym, meta, group_list_TN, file = exprSym_f)


