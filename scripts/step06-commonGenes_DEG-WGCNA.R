#To get the intersection of DEGs and WGCNA results with Veen diagram

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(VennDiagram)

lnames = load(file = file.path(Rdata_dir,'need_DEG_by_ImmuneScore.Rdata'))
lnames
lnames = load(file = paste(Rdata_dir,"WGCNA_03_ModuleBlueGenes_GS_MM.RData"))
lnames
commonGenes = intersect(rownames(need_DEG_by_ImmuneScore), rownames(ModuleBlue_GS_MM))
length(commonGenes)
commonGenes_logFC = need_DEG_by_ImmuneScore[commonGenes,]
commonGenes_GS_MM = ModuleBlue_GS_MM[commonGenes,
                                     c('GS.Immune_score','p.GSImmune_score','MMblue','p.MMblue')]
identical(rownames(commonGenes_logFC),rownames(commonGenes_GS_MM))
commonGenes_logFC_GS_MM = cbind(commonGenes_logFC,commonGenes_GS_MM)
save(commonGenes_logFC_GS_MM,
     file = paste(Rdata_dir,"commonGenes_DEG-WGCNA.RData"))

# venn plot
if(T){
  A = rownames(UP_DEG_by_ImmuneScore)
  B = rownames(DOWN_DEG_by_ImmuneScore)
  C = rownames(ModuleBlue_GS_MM)
  venn.plot <- venn.diagram(
    list("Up-DEGs"=A,"Down-DEGs"=B,"Blue module"=C),
    filename = NULL,
    lty = 1,
    lwd = 1,
    col = "black", 
    fill = c("red", "green", "blue"),
    alpha = 0.60,
    cex = 0.6, 
    reverse=F,
    rotation.degree = 0,
    cat.pos = 180, #position of each category name along the circle, 0 at 12 o'clock
    cat.dist = 0.12,#distance of each category name from the edge of the circle
    cat.col = "black",
    cat.cex = 0.6,# size for each category name
    cat.fontface ='plain'#fontface for each category name
  )
  dev.off()
  grid.draw(venn.plot)
  Fig_name=file.path(Figure_dir,'VennPlot_commonGenes.pdf')
  pdf(Fig_name,
      width = 2.7,
      height = 2.7)
  par(mai = c(0.25,0,0.25,0))
  grid.draw(venn.plot)
  dev.off()
}

#Indetify hub genes from commonGenes
#load hubGene from WGCNA
load(file = paste(Rdata_dir,"WGCNA_03_hubGene_ModuleBlue.RData"))
hubGenes = intersect(rownames(hubGene_ModuleBlue), commonGenes)
hubGenes_logFC_GS_MM = commonGenes_logFC_GS_MM[hubGenes,]
save(hubGenes_logFC_GS_MM,
     file = paste(Rdata_dir,"hubGenes_DEG-WGCNA.RData"))











