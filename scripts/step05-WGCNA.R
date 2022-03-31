#WGCNA to get modules related to immune score

options(stringsAsFactors = F)
Rdata_dir='../Rdata/'
Figure_dir='../figures/'
Table_dir='../tables/'
library(WGCNA)
enableWGCNAThreads()

#Step 01:Preparation of input data----
#Input data includes expression matrix (datExpr) and clinical data(datTraits)
if(T){
  workflow_type = "HTSeq - FPKM"
  exprSym_f=file.path(Rdata_dir,
                      paste(paste("TCGA",request_cancer,workflow_type,'exprSym_meta_groupT-N.RData',sep="-")))
  load(file = exprSym_f)
  dim(exprSym_FPKM)
  table(group_list_TN)
  exprSet=na.omit(exprSym_FPKM[,group_list_TN=='Tumor'])
  exprSet = log2(exprSet+1)
  dim(exprSet)
  
  #Selact genes with top 5000 variance for WGCNA
  datExpr = exprSet[order(apply(exprSet,1,var), decreasing = T)[1:5000],]
  datExpr = as.data.frame(t(datExpr)) 
  rownames(datExpr) = substr(rownames(datExpr),1,15)
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
}

#01:Test expression data
if(T){
  gsg = goodSamplesGenes(datExpr, verbose = 3);
  gsg$allOK
  if(!gsg$allOK){
    if(sum(!gsg$goodGenes)>0)
      table(gsg$goodGenes)
    printFlush(paste("removing genes:",paste(colnames(datExpr)[!gsg$goodGenes])))
    if(sum(!gsg$goodSamples)>0)
      table(gsg$goodSamples)
    printFlush(paste("removing samples:",paste(rownames(datExpr)[!gsg$goodSamples])))
    datExpr=datExpr[gsg$goodSamples,gsg$goodGenes]
    dim(datExpr) 
  }
  
  #Sample clustering to check outliers
  sampleTree = hclust(dist(datExpr), method = "average")
  Fig_name=file.path(Figure_dir, 'WGCNA-Sample clustering to detect outliers.pdf')
  pdf(Fig_name)
  par(mfrow=c(1,1))
  plot(sampleTree, main = "Sample clustering to detect outliers", 
       sub="", xlab="")
  dev.off()
  
  #Processing of outliers
  plot(sampleTree, main = "Sample clustering to detect outliers", 
       sub="", xlab="")
  cutHeight = 130
  abline(h=cutHeight,col="red")#Draw a line to show the culling criteria
  clust = cutreeStatic(sampleTree,cutHeight = cutHeight,minSize=10)
  table(clust)
  datExpr=datExpr[clust==1,]
  nGenes = ncol(datExpr)
  nGenes 
  nSamples = nrow(datExpr)
  nSamples 
}

#02:Processing sample information（Immune score）
if(T){
  #load estimate_score
  score_f=file.path(Rdata_dir,
                    paste(paste("TCGA",request_cancer,'tumor_estimateScore.RData',sep="-")))
  load(file = score_f)
  dim(estimate_score) 
  estimate_score = data.frame(estimate_score)
  
  traitRows=na.omit(match(rownames(datExpr),estimate_score$ID))
  datTraits=estimate_score[traitRows,]
  rownames(datTraits) = datTraits$ID 
  datTraits = datTraits[,c(2,3)]
  datExpr=datExpr[rownames(datTraits),] 
  identical(rownames(datTraits),rownames(datExpr))
  
  collectGarbage()
  #Recluster the samples上
  sampleTree2=hclust(dist(datExpr),method="average")
  #Translate clinical information into heatmaps: white-lowest,red-highest  
  traitColors=numbers2colors(datTraits,signed=FALSE)
  #cluster tree and heatmap
  plotDendroAndColors(sampleTree2,traitColors,
                      groupLabels=names(datTraits),
                      main="Sample dendrogram and trait heatmap")
}
#Save data
save(datExpr, datTraits, file = paste(Rdata_dir,"WGCNA_01_datInput.RData"))

#Step 02:Modules construction----
#screen soft threshold to obtain adjacency matrix, and then draw gene cluster tree to obtain Modules  
lnames=load(file.path(paste(Rdata_dir,"WGCNA_01_datInput.RData")))
lnames
#1.Screen soft threshold 
if(T){
  #soft-thresholding powers or beta
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  sft$powerEstimate 
  
  cex1 = 0.9
  dat = sft$fitIndices
  dat$R2 = -sign(dat[,3])*dat[,2]
  colnames(dat)[1]='Powers'
  library(ggplot2)
  #Association of “Scale Free Topology Fit index” and power
  if(T){
    p <- ggplot(dat, 
                aes(Powers, R2, label = Powers)) 
    p <- p + geom_text(size = 2.5,colour='red') + 
      geom_hline(aes(yintercept=cex1),lty=1,lwd=0.4,colour="red")+
      labs(title="Scale independence",
           x="Soft Threshold (power)",
           y="Scale Free Topology Model Fit,signed R^2") 
    p1_beta_Select <- p + theme(
      panel.background = element_rect(fill = 'white'), 
      plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5,face='bold'),
      plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
      plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
      axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
      axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
      axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
      axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
      legend.title  = element_text(color = 'black', size  = 6),
      legend.text   = element_text(color = 'black', size   = 6),
      axis.line = element_blank(), 
      panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
    )
    print(p1_beta_Select)
    save(p1_beta_Select,file = paste(Rdata_dir,"p1_WGCNA_beta_Select.RData"))
  }
  #Association of "mean connectivity" and power
  if(T){
    p <- ggplot(dat, 
                aes(Powers, mean.k., label = Powers)) 
    p <- p + geom_text(size = 2.5,colour='red') + 
      labs(title="Mean connectivity",
           x="Soft Threshold (power)",
           y="Mean Connectivity") 
    p2_beta_Select <- p + theme(
      panel.background = element_rect(fill = 'white'), 
      plot.title    = element_text(color = 'black', size   = 8.5, hjust = 0.5,face='bold'),
      plot.subtitle = element_text(color = 'black', size   = 6,hjust = 0.5),
      plot.caption  = element_text(color = 'black', size   = 6,face = 'italic', hjust = 1),
      axis.text.x   = element_text(color = 'black', size = 6, angle = 0),
      axis.text.y   = element_text(color = 'black', size = 6, angle = 0),
      axis.title.x  = element_text(color = 'black', size = 7.5, angle = 0),
      axis.title.y  = element_text(color = 'black', size = 7.5, angle = 90),
      legend.title  = element_text(color = 'black', size  = 6),
      legend.text   = element_text(color = 'black', size   = 6),
      axis.line = element_blank(), 
      panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA) 
    )
    print(p2_beta_Select)
    save(p2_beta_Select,file = paste(Rdata_dir,"p2_WGCNA_beta_Select.RData"))
  }
  #The two pictures above are drawn together
  p1_beta_Select <- p1_beta_Select + theme(
    plot.margin=unit(c(0.25, 0.2, 0.25, 0.25), 'cm'))#Adjust margins
  p2_beta_Select <- p2_beta_Select + theme(
    plot.margin=unit(c(0.25, 0.25, 0.25, 0.2), 'cm'))
  library(cowplot)
  plot_grid(p1_beta_Select, p2_beta_Select, 
            labels = NULL, 
            align = 'h')
  ggsave(filename = file.path(Figure_dir,'WGCNA-seletion_best_power_10-8cm.pdf'),
         plot = last_plot(), width = 10, height = 8, units = 'cm', dpi = 600)
  

  best_power = 3
  #Check whether the selected soft threshold is appropriate
  Connectivity = softConnectivity(datExpr,power=best_power)
  Fig_name=file.path(Figure_dir,paste0('WGCNA_bestPower-',best_power,'-scaleFreePlot.pdf'))
  pdf(Fig_name)
  par(mfrow=c(1,1))
  scaleFreePlot(Connectivity,main=paste("soft threshold,power=",best_power),truncated=T)
  dev.off()
}

#2.Network construction and module determination
if(T){
  #Use the following functions to generate the network in one step
  #If an error occurs when running the function, 
  #it may because of the confliction between WGCNA package and other functions  
  #You can try this:
  #       cor <- WGCNA::cor  #Firstly,Temporary assign function
  #       blockwiseModules； #Run blockwiseModules functions
  #       cor<-stats::cor
  
  #Basic process: 1)Calculate the adjacency between genes and then get the similarity between genes.
  # 2) derive the coefficient of intergene heterogeneity, and draw the systematic clustering tree.  
  # 3) Then set the minimum gene number of per module according to the standard of dynamic cut method.  
  # 4) calculate the characteristic vector values of each module.  
  # 5) cluster analysis to merge close modules.
  mergeHeight = 0.25
  net = blockwiseModules(datExpr, power = best_power, maxBlockSize = 8000,
                         TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = mergeHeight,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "WGCNA_gene_TOM",
                         verbose = 3)
  table(net$colors)
  
  #Draw the gene cluster tree to determine modules
  mergedColors = labels2colors(net$colors)# Convert labels to colors for plotting
  table(mergedColors)
  # Plot the dendrogram and the module colors underneath
  filename = file.path(Figure_dir,'WGCNA_gene_clusterDendrogram_and_module.pdf')
  pdf(filename) 
  plotDendroAndColors(net$dendrograms[[1]], 
                      mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

#Save data
if(T){
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  table(moduleColors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  save(MEs, moduleLabels, moduleColors, geneTree,
       file = paste(Rdata_dir,"WGCNA_02_networkConstruction.RData"))
}

#Step 03:Associate modules with sample information----
lnames=load(file.path(paste(Rdata_dir,"WGCNA_01_datInput.RData")))
lnames
lnames=load(file.path(paste(Rdata_dir,"WGCNA_02_networkConstruction.RData")))
lnames
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

#1.Average expression value of each module was calculated with moduleEigengenes  
if(T){
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  #To get gene significance
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
}
#2.Draw heatmap of the correlation between modules and sample information  
if(T){
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(", 
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  Fig_name=file.path(Figure_dir,'WGCNA_Module-trait_relationships.pdf')
  pdf(Fig_name,
      width = 3.5,#inch
      height = 4.5)
  par(mai = c(0.25,1.02,0.12,0.12))
  labeledHeatmap(Matrix = moduleTraitCor, 
                 xLabels = colnames(datTraits), 
                 yLabels = colnames(MEs), 
                 ySymbols = colnames(MEs),
                 xLabelsAngle = 0,
                 xLabelsAdj = 0.5, 
                 cex.lab = 0.75, 
                 yColorWidth=0.01, 
                 xColorWidth = 0.01,
                 colorLabels = FALSE, 
                 colors = blueWhiteRed(50), 
                 textMatrix = textMatrix, 
                 setStdMargins = FALSE, 
                 cex.text = 0.7, 
                 zlim = c(-1,1) 
                 )
  dev.off()
}
#3.To calculate ModuleMembership
if(T){
  modNames=substring(colnames(MEs),3)
  geneModuleMembership=as.data.frame(cor(datExpr,MEs,use="p"))
  MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  colnames(geneModuleMembership)=paste("MM",modNames,sep="")
  colnames(MMPvalue)=paste("p.MM",modNames,sep="")
}

#4.Calculate the correlation coefficient between gene and the clinical feature of interest（GS，gene significance）
if(T){
  geneTraitSignificance=as.data.frame(cor(datExpr,datTraits,use="p"))
  GSPvalue=as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  colnames(geneTraitSignificance)=paste("GS.",colnames(datTraits),sep="")
  colnames(GSPvalue)=paste("p.GS",colnames(datTraits),sep="")
}

#Extract genes contained in the module of interest
if(T){
  identical(rownames(geneTraitSignificance),rownames(geneModuleMembership))
  identical(rownames(GSPvalue),rownames(MMPvalue))
  identical(rownames(geneTraitSignificance),rownames(GSPvalue))
  gene_GS_MM = cbind(geneTraitSignificance,GSPvalue,geneModuleMembership,MMPvalue)
  
  module="blue"
  column=match(module,modNames)
  moduleGenes=(moduleColors==module)
  table(moduleGenes)
  ModuleBlue_GS_MM = gene_GS_MM[moduleGenes,] 
  save(ModuleBlue_GS_MM,
       file = paste(Rdata_dir,"WGCNA_03_ModuleBlueGenes_GS_MM.RData"))
}

#Screen hub genes based on MM and GS(|MM|>0.9,|GS|>0.8)
if(T){
  t_MM = 0.9
  t_GS = 0.8
  t_pvalue = 0.01
  identical(rownames(geneTraitSignificance),rownames(geneModuleMembership))
  identical(rownames(GSPvalue),rownames(MMPvalue))
  identical(rownames(geneTraitSignificance),rownames(GSPvalue))
  gene_GS_MM = cbind(geneTraitSignificance,GSPvalue,geneModuleMembership,MMPvalue)
  
  module="blue"
  column=match(module,modNames)
  moduleGenes=(moduleColors==module)
  table(moduleGenes)
  ModuleBlue_GS_MM = gene_GS_MM[moduleGenes,] 
  save(ModuleBlue_GS_MM, file = paste(Rdata_dir,"WGCNA_03_ModuleBlueGenes_GS_MM.RData"))
  hubGene_ModuleBlue = ModuleBlue_GS_MM[abs(ModuleBlue_GS_MM$GS.Immune_score)>t_GS
                                        &(ModuleBlue_GS_MM$p.GSImmune_score<t_pvalue)
                                        &(abs(ModuleBlue_GS_MM$MMblue)>t_MM)
                                        &(ModuleBlue_GS_MM$p.MMblue<t_pvalue),]
  write.csv(hubGene_ModuleBlue[,c('MMblue','p.MMblue','GS.Immune_score','p.GSImmune_score')],
            file = file.path(Table_dir,'WGCNA_hubGenes_MM-GS.csv'),
            col.names = T, row.names = T)
  save(hubGene_ModuleBlue,
       file = paste(Rdata_dir,"WGCNA_03_hubGene_ModuleBlue.RData"))
  
  #scatter plot of correlation between GS and MM
  if(T){
    dat = ModuleBlue_GS_MM[,c('MMblue','GS.Immune_score')]
    dat$hub = ifelse(rownames(dat) %in% rownames(hubGene_ModuleBlue),'T','F')
    table(dat$hub) 
    library(ggpubr)
    p <- ggplot(data=dat,
                aes(x=MMblue, y=GS.Immune_score, color=hub)) +
      geom_point(alpha=0.8, size=0.4) + 
      scale_colour_manual(values = c('grey','blue'))+ 
      xlim(-0.5,1.0)+ 
      ylim(-0.5,1.0)+
      xlab("MM in blue module") +
      ylab("GS for immune score") +
      geom_hline(yintercept = t_GS,lty=2,lwd=0.4,alpha=0.8) +
      geom_vline(xintercept = t_MM,lty=2,lwd=0.4,alpha=0.8)+
      theme(
        panel.background = element_rect(fill = 'white'), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        plot.title    = element_text(color = 'black', size   = 12, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 7.5,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 7.5,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 7.5, angle = 0),
        axis.text.y   = element_text(color = 'black', size = 7.5, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 10, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 10, angle = 90),
        legend.title  = element_blank(),
        legend.text   = element_text(color = 'black', size   = 7.5),
        legend.position="none",
        axis.line = element_blank(), 
        panel.border = element_rect(linetype = 'solid', size = 0.6,fill = NA)
      )
    print(p)
    ggsave(filename = file.path(Figure_dir,'cor_GS_MM.pdf'),
           plot = last_plot(), width = 7.9, height = 7.9, units = 'cm', dpi = 300)
  }
}

#5.Further test whether the selected Module is indeed related to the clinical features by scatter plots  
if(T){
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                     abs(geneTraitSignificance[moduleGenes,1]),
                     xlab=paste("Module Membership in",module,"module"),
                     ylab = "Gene significance for stromal score",
                     main = paste("Module Membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                     abs(geneTraitSignificance[moduleGenes,2]),
                     xlab=paste("Module Membership in",module,"module"),
                     ylab = "Gene significance for immune score",
                     main = paste("Module Membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
}

#Step 04:Result visualization  ----
lnames=load(file.path(paste(Rdata_dir,"WGCNA_01_datInput.RData")))
lnames
lnames=load(file.path(paste(Rdata_dir,"WGCNA_02_networkConstruction.RData")))
lnames
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

#1.TOM network
if(T){
  # Calculate topological overlap anew: this could be done more efficiently by saving the TOM
  # calculated during module detection, but let us do it again here.
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = best_power)
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM^7
  dim(plotTOM)
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA

  #TOM network with 400 genes
  if(T){
    nSelect = 1000
    # For reproducibility, we set the random seed
    set.seed(10)
    select = sample(nGenes, size = nSelect)
    selectTOM = dissTOM[select, select]
    # There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
    selectTree = hclust(as.dist(selectTOM), method = "average")
    selectColors = moduleColors[select]
    # Open a graphical window
    # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
    # the color palette; setting the diagonal to NA also improves the clarity of the plot
    plotDiss = selectTOM^7;
    diag(plotDiss) = NA;
    Fig_name=file.path(Figure_dir,'WGCNA_TOM_plot.pdf')
    pdf(Fig_name,
        width = 3.5,
        height = 3.5)
    TOMplot(plotDiss, selectTree, 
            selectColors, 
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
    dev.off()
  }
}

# 2.draw the relationship between modules and traits
if(T){
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  CB = as.data.frame(design[,2]);
  names(CB) = "CB"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, CB))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  par(cex = 0.9)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  # Plot the dendrogram
  sizeGrWindow(6,6);
  par(cex = 1.0)
  plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90)
}

#Export gene co-expression network of the specified modules,
#which can be further imported into Cytoscape to draw the network diagram  
if(T){
  lnames=load(file.path(paste(Rdata_dir,"WGCNA_01_datInput.RData")))
  lnames
  lnames=load(file.path(paste(Rdata_dir,"WGCNA_02_networkConstruction.RData")))
  lnames
  lnames=load(file = paste(Rdata_dir,"WGCNA_03_hubGene_ModuleBlue.RData"))
  lnames
  
  TOM = TOMsimilarityFromExpr(datExpr, power = best_power)
  modules = c("blue")
  probes = names(datExpr)
  inModule = is.finite(match(moduleColors, modules))
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  rownames(modTOM) = modProbes
  colnames(modTOM) = modProbes
  hubGeneTOM = modTOM[rownames(hubGene_ModuleBlue),rownames(hubGene_ModuleBlue)]
  
  cyt = exportNetworkToCytoscape(
    hubGeneTOM, 
    edgeFile = paste(Table_dir,"CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
    nodeFile = paste(Table_dir,"CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.02,#adjacency threshold for including edges in the output
    nodeNames = rownames(hubGeneTOM)
  )
  #Two files are obtained finally to import into cytoscape:
  #CytoscapeInput-edges-brown-red.txt
  #CytoscapeInput-nodes-brown-red.txt
}



