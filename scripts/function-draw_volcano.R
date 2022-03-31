#Draw volcano based on DEGs

draw_volcano <- function(DEG_result,pvalue_t,logFC_t,method,figure_dir){
  # Note: DEG_result is the result of DEG,the format is dataframe,rownames--gene,2 columns--log2FoldChange,pvalue.
  # pvalue_t,logFC_t is the threshold for significance.
  # method, such as 'DEseq2','edgeR','limma'.
  # figure_dir is the path to save figures,for example, figure_dir="../figures/"
  
  library(ggplot2)
  DEG_result$change= (ifelse(DEG_result$pvalue < pvalue_t & abs(DEG_result$log2FoldChange) >logFC_t,
                             ifelse(DEG_result$log2FoldChange > logFC_t,'Up','Down'),'Unchanged'))
  DEG_result = DEG_result[DEG_result$pvalue>0,]
  table(DEG_result$change)
  p_volcano = ggplot(data=DEG_result,
             aes(x=log2FoldChange, y=-log10(pvalue), color=change)) +
    geom_point(alpha=0.6, size=0.6) + #alpha is transparency of dots
    scale_colour_manual(values = c('blue','grey','red'))+ 
    xlim(-3,5)+ 
    xlab("log2(FoldChange)") +
    ylab("-log10(FDR)") +
    geom_hline(yintercept = -log10(pvalue_t),lty=2,lwd=0.4,alpha=0.8) +
    geom_vline(xintercept = c(logFC_t,(0-logFC_t)),lty=2,lwd=0.4,alpha=0.8)+
    theme(
      panel.background = element_rect(fill = 'white'), 
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent"),
      plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5),
      plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
      plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
      axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
      axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
      axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
      axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
      legend.title  = element_blank(),
      legend.text   = element_text(color = 'black', size   = 6),
      #legend.position=c(0.14,0.87),
      legend.position="top",
      axis.line = element_blank(), 
      #axis.line.y = element_line(color = 'black', linetype = 'solid'), 
      #axis.line.x = element_line (color = 'black',linetype = 'solid'), 
      panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
    )
  print(p_volcano)
  save(p_volcano, file = file.path(Rdata_dir,'p_volcano.Rdata'))
  ggsave(filename = file.path(figure_dir,paste0(method,'_DEG_volcano_7-9cm.pdf')),
         plot = last_plot(), width = 7, height = 9, units = 'cm', dpi = 600)
}
