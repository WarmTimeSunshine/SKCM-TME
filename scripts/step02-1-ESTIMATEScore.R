#Read downloaded estimate score of TCGA data

#rm(list = ls()) 
options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'

if(T){
  estimate_score = read.table(file = file.path(Table_dir, paste("ESTIMATE",request_cancer,'fromWeb.txt',sep="_")),
                              header=TRUE)
  dim(estimate_score) #471 tumor samples
  
  estimate_score = data.frame(estimate_score)
  table(table(estimate_score$ID)) 
  dim(estimate_score)
  
  #Save data----
  score_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,'tumor_estimateScore.RData',sep="-")))
  save(estimate_score, file = score_f)
}



