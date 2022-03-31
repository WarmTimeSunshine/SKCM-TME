#\GO,KEGG enrichment 

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
#load genes for enrichment analysis---commonGenes
load(file = paste(Rdata_dir,"commonGenes_DEG-WGCNA.RData"))

#Transform genes names from symbol to ENTREZID
if(T){
  geneName <- bitr(rownames(commonGenes_logFC_GS_MM), fromType = "SYMBOL",
                   toType = c( "ENTREZID"),
                   OrgDb = org.Hs.eg.db)
}

options(digits=3) 
#GO enrichment 
if(T){
  if(T){
    ego_all <- enrichGO(gene = geneName$ENTREZID, 
                        OrgDb = "org.Hs.eg.db", 
                        ont="all",
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = 0.01,
                        qvalueCutoff =0.05,
                        readable      = TRUE)
    ego_all=DOSE::setReadable(ego_all, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    ego_all_dt=ego_all@result
    save(ego_all_dt,file = file.path(Rdata_dir,'_GO_ego_all_dt_results.RData'))
    write.csv(ego_all_dt,file = file.path(Table_dir,"enrich_GO_commonGenes.csv"),
              col.names = TRUE, row.names = TRUE)
    
    #Draw bar chart
    if(F){
      enrich_BP <- ego_all_dt[ego_all_dt$ONTOLOGY=='BP',] 
      enrich_BP <- enrich_BP[order(enrich_BP$p.adjust),]
      enrich_BP <- enrich_BP[1:10,]  #top10 pathways
      enrich_BP$Description[1:3] = c("immune-activating receptor signaling pathway",
                                     "immune-activating signal transduction",
                                     "adaptive immune response")
      count <- as.numeric(unlist(strsplit(enrich_BP$GeneRatio,"/655",fixed=T))) 
      enrich_BP <- data.frame(enrich_BP[,1:3],count,enrich_BP$p.adjust)
      colnames(enrich_BP)[5] <- c("p.adj")
      
      enrich_CC <- ego_all_dt[ego_all_dt$ONTOLOGY=='CC',] 
      enrich_CC <- enrich_CC[order(enrich_CC$p.adjust),] 
      enrich_CC <- enrich_CC[1:11,]  
      enrich_CC <- enrich_CC[-9,]  
      enrich_CC$Description[9] = "lumenal side of ER membrane"
      count <- as.numeric(unlist(strsplit(enrich_CC$GeneRatio,"/678",fixed=T))) 
      enrich_CC <- data.frame(enrich_CC[,1:3],count,enrich_CC$p.adjust)
      colnames(enrich_CC)[5] <- c("p.adj")
      
      enrich_MF <- ego_all_dt[ego_all_dt$ONTOLOGY=='MF',] 
      enrich_MF <- enrich_MF[order(enrich_MF$p.adjust),] 
      enrich_MF <- enrich_MF[1:10,] 
      count <- as.numeric(unlist(strsplit(enrich_MF$GeneRatio,"/626",fixed=T))) 
      enrich_MF <- data.frame(enrich_MF[,1:3],count,enrich_MF$p.adjust)
      colnames(enrich_MF)[5] <- c("p.adj")
      
      enrich3 = rbind(enrich_BP,enrich_CC,enrich_MF)
      library(Hmisc)
      enrich3$Description = capitalize(enrich3$Description)
      p <- ggplot(data=enrich3,
                  aes(x=factor(Description,levels = rev(enrich3$Description)),
                      y=count,fill=p.adj)) 
      p1 <- p + geom_bar(stat="identity",position = "dodge") + 
        coord_flip() + 
        scale_fill_gradient(low="red",high="blue") +   
        labs(x=NULL,y="Gene Counts") +
        facet_grid(enrich3$ONTOLOGY~.,scales = "free",space = "free")
      p_boxplot_GO <- p1 + theme(
        panel.background=element_rect(fill='transparent',color='gray'),
        plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
        axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
        axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
        axis.line = element_blank(), 
        panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
      )
      print(p_boxplot_GO)
      save(p_boxplot_GO, file = file.path(Rdata_dir,'p_boxplot_GO.Rdata'))
      ggsave(filename = file.path(Figure_dir,'GO_all_boxplot.pdf'),
             plot = last_plot(), width = 18, height = 14, units = 'cm', dpi = 300)
    }
  }
  
  save(ego_BP_dt,ego_MF_dt,ego_CC_dt,
       file = file.path(Rdata_dir,'_GO_ego_all_dt_results.RData'))
}

## KEGG 
if(T){
  kegg_result <- enrichKEGG(gene  = geneName$ENTREZID,
                            organism     = 'hsa',
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.01,
                            minGSSize = 10,
                            maxGSSize = 500)
  kegg_dt <- as.data.frame(kegg_result)
  save(kegg_dt,file = file.path(Rdata_dir,'_KEGG_enrich_results.RData'))
  write.csv(kegg_dt,file = file.path(Table_dir,"enrich_KEGG_commonGenes.csv"),
            col.names = TRUE, row.names = TRUE)
  
  #dotplot
  dotplot(kegg_result,showCategory = 20)
  ggsave(file.path(Figure_dir,'KEGG_ego_all_dt_dotplot.pdf'))
  
  #barplot
  barplot(kegg_result,showCategory = 10) +
    ggsave(file.path(Figure_dir,'KEGG_ego_all_dt_barplot.pdf'))
}

#Draw the results of GO and KEGG together
if(T){
  if(T){
    enrich_BP <- ego_all_dt[ego_all_dt$ONTOLOGY=='BP',] 
    enrich_BP <- enrich_BP[order(enrich_BP$p.adjust),] 
    enrich_BP <- enrich_BP[1:10,]  
    enrich_BP$Description[c(1,3,10)] = c("immune receptor mediated adaptive immune response",
                                         "immune receptor signaling pathway",
                                         "immunoglobulin mediated humoral immune response")
    count <- as.numeric(unlist(strsplit(enrich_BP$GeneRatio,"/578",fixed=T))) 
    enrich_BP <- data.frame(enrich_BP[,1:3],count,enrich_BP$p.adjust)
    colnames(enrich_BP)[5] <- 'p.adjust'
    
    enrich_CC <- ego_all_dt[ego_all_dt$ONTOLOGY=='CC',] 
    enrich_CC <- enrich_CC[order(enrich_CC$p.adjust),] 
    enrich_CC <- enrich_CC[1:11,]  
    enrich_CC <- enrich_CC[-6,]  
    count <- as.numeric(unlist(strsplit(enrich_CC$GeneRatio,"/601",fixed=T))) 
    enrich_CC <- data.frame(enrich_CC[,1:3],count,enrich_CC$p.adjust)
    colnames(enrich_CC)[5] <- 'p.adjust'
    
    enrich_MF <- ego_all_dt[ego_all_dt$ONTOLOGY=='MF',] 
    enrich_MF <- enrich_MF[order(enrich_MF$p.adjust),] 
    enrich_MF <- enrich_MF[1:10,]  
    count <- as.numeric(unlist(strsplit(enrich_MF$GeneRatio,"/575",fixed=T))) 
    enrich_MF <- data.frame(enrich_MF[,1:3],count,enrich_MF$p.adjust)
    colnames(enrich_MF)[5] <- 'p.adjust'
    
    enrich_GO = rbind(enrich_BP,enrich_CC,enrich_MF)
    library(Hmisc)
    enrich_GO$ONTOLOGY = paste0('GO_',enrich_GO$ONTOLOGY)
    enrich_GO$Description = capitalize(enrich_GO$Description)
    
    enrich_KEGG <- kegg_dt[order(kegg_dt$p.adjust),]
    enrich_KEGG <- enrich_KEGG[1:10,]
    enrich_KEGG <- cbind(rep('KEGG',nrow(enrich_KEGG)),enrich_KEGG)
    colnames(enrich_KEGG)[1] <- 'ONTOLOGY'
    count <- as.numeric(unlist(strsplit(enrich_KEGG$GeneRatio,"/323",fixed=T))) 
    enrich_KEGG <- data.frame(enrich_KEGG[,1:3],count,enrich_KEGG$p.adjust)
    colnames(enrich_KEGG)[5] <- 'p.adjust'
    
    enrich3 = rbind(enrich_GO,enrich_KEGG)
  }
  
  #bar plot
  if(T){
    p <- ggplot(data=enrich3,
                aes(x=factor(Description,levels = rev(enrich3$Description)),
                    y=count,fill=p.adjust)) 
    p1 <- p + geom_bar(stat="identity",position = "dodge") + 
      coord_flip() +  
      scale_fill_gradient(low="red",high="blue") +   
      labs(x=NULL,y="Gene Counts") +
      facet_grid(enrich3$ONTOLOGY~.,scales = "free",space = "free")
    p_boxplot_GO_KEGG <- p1 + theme(
      panel.background=element_rect(fill='transparent',color='gray'),
      plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
      axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
      axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
      axis.line = element_blank(), 
      panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
    )
    print(p_boxplot_GO_KEGG)
    save(p_boxplot_GO_KEGG, file = file.path(Rdata_dir,'p_boxplot_GO_KEGG.Rdata'))
    ggsave(filename = file.path(Figure_dir,'boxplot_GO_KEGG.pdf'),
           plot = last_plot(), width = 18, height = 14, units = 'cm', dpi = 300)
  }
}



