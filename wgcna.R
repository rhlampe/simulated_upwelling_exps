#
# WGCNA boilerplate for analyzing UCBC responses
# 
# This code is largely based on the WGCNA tutorials:
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA

# Load library
library(WGCNA)

#
# Read in normalized/filtered count data
#
vsd_counts <- read.csv('vsd_counts_filtered.csv', row.names = 1)
vsd_counts <- as.data.frame(t(vsd_counts)) # goodSampleGenes requires this format

#
# Load 'trait' data. This should correspond to 1 or 0 for UCBC stage
#
datTraits <- read.csv('sampleData.csv', row.names = 1)

#
# Clean up count data
#
datExpr <- vsd_counts

# Run goodSampleGenes()
gsg = goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
	printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
	if (sum(!gsg$goodSamples)>0)
		printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
	datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

#
# cluster samples based on Euclidean distance and look at dendrogram for outliers
#
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree)

#
# Select soft power threshold
#
powers = c(c(1:10), seq(from = 10, to = 30, by = 1)) # choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed") # call network topology analysis function

sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# Set soft power threshold based on above plots
softPower <- NA

#
# Construct Networks
#
enableWGCNAThreads()
gc()
adjacency = adjacency(datExpr, power = softPower, type = "signed")

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType = "signed") # specify network type
dissTOM = 1-TOM

#
# Generate a clustered gene tree
#
library(flashClust)

geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")

#plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
#set a threhold for merging modules. Here we are not merging so MEDissThres=0.0
MEDissThres = 0.0
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")

#
# Correlate traits
#

#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))

#display the corelation values with a heatmap plot
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))


#
# Gene module membership and significance
#
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


# Example for finding significant genes for time point 1. Modify or repeat for other time points
tp1 <- as.data.frame(datTraits$Tp_1)

geneTraitSignificance = as.data.frame(cor(datExpr, tp1, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- "GS.tp1"
names(GSPvalue) <- "p.GS.tp1"

# Retrieve genes with trait significance and module membership
tp1_mes <- gsub("ME", "", tp_1)
tp1_out <- NA
for (me in tp1_mes) {
  tp1_genes <- names(datExpr)[moduleColors==me]
  mm <- geneModuleMembership[tp1_genes, paste("MM", me, sep="")]
  mmp <- MMPvalue[tp1_genes,paste("p.MM", me, sep="")]
  mmp.adj <- p.adjust(mmp, method="BH")
  gs <- geneTraitSignificance[tp1_genes,1]
  gsp <- GSPvalue[tp1_genes,1]
  gsp.adj <- p.adjust(gsp, method="BH")
  
  if (is.na(tp1_out)) {
    tp1_out <- data.frame(tp1_genes, mm, mmp.adj, gs, gsp.adj)
    tp1_out <- tp1_out[tp1_out$mmp.adj < 0.05,]
    tp1_out <- tp1_out[tp1_out$gsp.adj < 0.05,]
  } else {
    temp_out <- data.frame(tp1_genes, mm, mmp.adj, gs, gsp.adj)
    temp_out <- temp_out[temp_out$mmp.adj < 0.05,]
    temp_out <- temp_out[temp_out$gsp.adj < 0.05,]
    
    tp1_out <- rbind(tp1_out, temp_out)
  }
}

write.csv(tp1_out, 'tp1_significant_genes.csv')
