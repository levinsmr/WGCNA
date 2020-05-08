#=====================================================================================
#
#  Import IP Data
#
#=====================================================================================
rm(list=ls())

# If necessary, change the path below to the directory where the data files are stored. 
workingDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/Data Sheets";
setwd(workingDir); 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Read in the data (Here I am reading Data that has been Mean Filterd)
Data = read.csv("GeneTableIP.csv");
# Data = read.csv("GeneTableIP_Mean.csv");
# Data = read.csv("GeneTableIP_Variance.csv");

# Set Up Gene Table
GeneCounts = as.data.frame(t(Data[, -c(1)]));
names(GeneCounts) = Data$Genes;

#=====================================================================================
#
#  Load WGCNA
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
enableWGCNAThreads(6)

#=====================================================================================
#
# Function checks data for missing entries and zero-variance genes, and returns a list of samples and genes that pass criteria maximum number of missing values. If necessary, the filtering is iterated.
#
#=====================================================================================

dim(GeneCounts);
gsg = goodSamplesGenes(GeneCounts, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(GeneCounts)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(GeneCounts)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  GeneCounts = GeneCounts[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Choose Soft Thresholding Power
#
#=====================================================================================
# Log TransformNormalization
# GeneCounts = log2(GeneCounts + 1);

# Choose a set of soft-thresholding powers
powers = c(c(1:30))

# Call the network topology analysis function
sft = pickSoftThreshold(GeneCounts, powerVector = powers, verbose = 3)

# Figure Directory
FigDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/WGCNA Output/Figures";
setwd(FigDir); 

# Plot the results:
sizeGrWindow(9, 5)
cex1 = 0.75;
png("Soft Power Threshold IP.png");
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
dev.off()

# Plot the results:
sizeGrWindow(9, 5)
cex1 = 0.75;
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#=====================================================================================
#
#  manual network construction and module detection
#
#=====================================================================================

adjacency_signed = adjacency(GeneCounts, power= 3, type="signed", corFnc="bicor");
TOM_signed = TOMsimilarity(adjacency_signed, TOMType="signed", verbose=3);
dissTOM_signed = 1-TOM_signed;
geneTree = hclust(as.dist(dissTOM_signed), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Morphine RNA-Seq IP", labels = FALSE, hang = 0.04);
dynamicMods_hybrid = cutreeDynamic(dendro = geneTree, distM = dissTOM_signed, method="hybrid", deepSplit = 3, pamRespectsDendro = TRUE, minClusterSize = 50);
moduleColor = labels2colors(dynamicMods_hybrid)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, moduleColor, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Morphine Microglia RNA-Seq")

#=====================================================================================
#
#  # Calculate eigengenes & Cluster Modules
#
#=====================================================================================

# Calculate eigengenes
MEList = moduleEigengenes(GeneCounts, moduleColor)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Plot the cut line into the dendrogram
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(GeneCounts, moduleColor, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(18, 9)
png("Merged gene Dendrogram IP.png")
plotDendroAndColors(geneTree, cbind(moduleColor, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
#
#  Save Data
#
#=====================================================================================

# Figure Directory
DataDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/WGCNA Output";
setwd(DataDir); 

# Save dissTOM
write.csv(dissTOM_signed, file = "dissTOM IP.csv",  row.names=FALSE)

# Save Module Eigengenes
# write.csv(MEs, file = "MEs IP.csv",  row.names=FALSE)

# Save Module Hub Genes
hubs = chooseTopHubInEachModule(GeneCounts, moduleColor)
# write.csv(hubs, file = "Hubs IP.csv",  row.names=FALSE)

# Save Module Membership & Gene Count
geneInfo = data.frame(names(GeneCounts),moduleColor,t(GeneCounts));
write.csv(geneInfo, file = "Gene Modules IP.csv")

# Rename to moduleColors
# moduleColor = mergedColors;

# Save Module Membership & Gene Count
# geneInfo = data.frame(names(GeneCounts),mergedColors,t(GeneCounts));
# write.csv(geneInfo, file = "Gene Modules Merged IP DeSeq_Minus.csv")

# Save Module Hub Genes
# hubs = chooseTopHubInEachModule(GeneCounts, mergedColors)
# write.csv(hubs, file = "Hubs Merged IP DeSeq_Minus.csv",  row.names=FALSE)

# Save Module Eigengenes Merged
# write.csv(mergedMEs, file = "MEs Merged IP DeSeq_Minus.csv",  row.names=FALSE)

#=====================================================================================
#
#  Import IP Adjusted for Noise Data
#
#=====================================================================================
rm(list=ls())

# If necessary, change the path below to the directory where the data files are stored. 
workingDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/Data Sheets";
setwd(workingDir); 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Read in the data (Here I am reading Data that has been Mean Filterd)
Data = read.csv("GeneTableIP_Minus.csv");
# Data = read.csv("GeneTableIP_Mean.csv");
# Data = read.csv("GeneTableIP_Variance.csv");

# Set Up Gene Table
GeneCounts = as.data.frame(t(Data[, -c(1)]));
names(GeneCounts) = Data$Genes;

#=====================================================================================
#
#  Load WGCNA
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
enableWGCNAThreads(6)

#=====================================================================================
#
# Function checks data for missing entries and zero-variance genes, and returns a list of samples and genes that pass criteria maximum number of missing values. If necessary, the filtering is iterated.
#
#=====================================================================================

dim(GeneCounts);
gsg = goodSamplesGenes(GeneCounts, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(GeneCounts)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(GeneCounts)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  GeneCounts = GeneCounts[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Choose Soft Thresholding Power
#
#=====================================================================================
# Log TransformNormalization
# GeneCounts = log2(GeneCounts + 1);

# Choose a set of soft-thresholding powers
powers = c(c(1:30))

# Call the network topology analysis function
sft = pickSoftThreshold(GeneCounts, powerVector = powers, verbose = 3)

# Figure Directory
FigDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/WGCNA Output/Figures";
setwd(FigDir); 

# Plot the results:
sizeGrWindow(9, 5)
cex1 = 0.75;
png("Soft Power Threshold IP_Minus.png");
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
dev.off()

# Plot the results:
sizeGrWindow(9, 5)
cex1 = 0.75;
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#=====================================================================================
#
#  manual network construction and module detection
#
#=====================================================================================

adjacency_signed = adjacency(GeneCounts, power= 3, type="signed", corFnc="bicor");
TOM_signed = TOMsimilarity(adjacency_signed, TOMType="signed", verbose=3);
dissTOM_signed = 1-TOM_signed;
geneTree = hclust(as.dist(dissTOM_signed), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Morphine RNA-Seq IP Minus", labels = FALSE, hang = 0.04);
dynamicMods_hybrid = cutreeDynamic(dendro = geneTree, distM = dissTOM_signed, method="hybrid", deepSplit = 3, pamRespectsDendro = TRUE, minClusterSize = 50);
moduleColor = labels2colors(dynamicMods_hybrid)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, moduleColor, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Morphine Microglia RNA-Seq")

#=====================================================================================
#
#  # Calculate eigengenes & Cluster Modules
#
#=====================================================================================

# Calculate eigengenes
MEList = moduleEigengenes(GeneCounts, moduleColor)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Plot the cut line into the dendrogram
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(GeneCounts, moduleColor, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(18, 9)
png("Merged gene Dendrogram IP Minus.png")
plotDendroAndColors(geneTree, cbind(moduleColor, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
#
#  Save Data
#
#=====================================================================================

# Figure Directory
DataDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/WGCNA Output";
setwd(DataDir); 

# Save dissTOM
write.csv(dissTOM_signed, file = "dissTOM IP Minus.csv",  row.names=FALSE)

# Save Module Eigengenes
# write.csv(MEs, file = "MEs IP.csv",  row.names=FALSE)

# Save Module Hub Genes
hubs = chooseTopHubInEachModule(GeneCounts, moduleColor)
# write.csv(hubs, file = "Hubs IP.csv",  row.names=FALSE)

# Save Module Membership & Gene Count
geneInfo = data.frame(names(GeneCounts),moduleColor,t(GeneCounts));
write.csv(geneInfo, file = "Gene Modules IP Minus.csv")

# Rename to moduleColors
# moduleColor = mergedColors;

# Save Module Membership & Gene Count
# geneInfo = data.frame(names(GeneCounts),mergedColors,t(GeneCounts));
# write.csv(geneInfo, file = "Gene Modules Merged IP DeSeq_Minus.csv")

# Save Module Hub Genes
# hubs = chooseTopHubInEachModule(GeneCounts, mergedColors)
# write.csv(hubs, file = "Hubs Merged IP DeSeq_Minus.csv",  row.names=FALSE)

# Save Module Eigengenes Merged
# write.csv(mergedMEs, file = "MEs Merged IP DeSeq_Minus.csv",  row.names=FALSE)

#=====================================================================================
#
#  Import Input Data
#
#=====================================================================================
rm(list=ls())

# If necessary, change the path below to the directory where the data files are stored. 
workingDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/Data Sheets";
setwd(workingDir); 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Read in the data (Here I am reading Data that has been Mean Filterd)
Data = read.csv("GeneTableIN.csv");
# Data = read.csv("GeneTableIP_Mean.csv");
# Data = read.csv("GeneTableIP_Variance.csv");

# Set Up Gene Table
GeneCounts = as.data.frame(t(Data[, -c(1)]));
names(GeneCounts) = Data$Genes;

#=====================================================================================
#
#  Load WGCNA
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
enableWGCNAThreads(6)

#=====================================================================================
#
# Function checks data for missing entries and zero-variance genes, and returns a list of samples and genes that pass criteria maximum number of missing values. If necessary, the filtering is iterated.
#
#=====================================================================================

dim(GeneCounts);
gsg = goodSamplesGenes(GeneCounts, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(GeneCounts)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(GeneCounts)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  GeneCounts = GeneCounts[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Choose Soft Thresholding Power
#
#=====================================================================================
# Log TransformNormalization
# GeneCounts = log2(GeneCounts + 1);

# Choose a set of soft-thresholding powers
powers = c(c(1:30))

# Call the network topology analysis function
sft = pickSoftThreshold(GeneCounts, powerVector = powers, verbose = 3)

# Figure Directory
FigDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/WGCNA Output/Figures";
setwd(FigDir); 

# Plot the results:
sizeGrWindow(9, 5)
cex1 = 0.75;
png("Soft Power Threshold IN.png");
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
dev.off()

# Plot the results:
sizeGrWindow(9, 5)
cex1 = 0.75;
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#=====================================================================================
#
#  manual network construction and module detection
#
#=====================================================================================

adjacency_signed = adjacency(GeneCounts, power= 5, type="signed", corFnc="bicor");
TOM_signed = TOMsimilarity(adjacency_signed, TOMType="signed", verbose=3);
dissTOM_signed = 1-TOM_signed;
geneTree = hclust(as.dist(dissTOM_signed), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Morphine RNA-Seq IN", labels = FALSE, hang = 0.04);
dynamicMods_hybrid = cutreeDynamic(dendro = geneTree, distM = dissTOM_signed, method="hybrid", deepSplit = 3, pamRespectsDendro = TRUE, minClusterSize = 50);
moduleColor = labels2colors(dynamicMods_hybrid)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, moduleColor, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Morphine Microglia RNA-Seq")

#=====================================================================================
#
#  # Calculate eigengenes & Cluster Modules
#
#=====================================================================================

# Calculate eigengenes
MEList = moduleEigengenes(GeneCounts, moduleColor)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Plot the cut line into the dendrogram
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(GeneCounts, moduleColor, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(18, 9)
png("Merged gene Dendrogram IN.png")
plotDendroAndColors(geneTree, cbind(moduleColor, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
#
#  Save Data
#
#=====================================================================================

# Figure Directory
DataDir = "Z:/Neumaier Lab/Morphine Grant/RNA Seq 2/WGCNA for Publication/WGCNA Output";
setwd(DataDir); 

# Save dissTOM
write.csv(dissTOM_signed, file = "dissTOM IN.csv",  row.names=FALSE)

# Save Module Eigengenes
# write.csv(MEs, file = "MEs IP.csv",  row.names=FALSE)

# Save Module Hub Genes
hubs = chooseTopHubInEachModule(GeneCounts, moduleColor)
# write.csv(hubs, file = "Hubs IP.csv",  row.names=FALSE)

# Save Module Membership & Gene Count
geneInfo = data.frame(names(GeneCounts),moduleColor,t(GeneCounts));
write.csv(geneInfo, file = "Gene Modules IN.csv")

# Rename to moduleColors
# moduleColor = mergedColors;

# Save Module Membership & Gene Count
# geneInfo = data.frame(names(GeneCounts),mergedColors,t(GeneCounts));
# write.csv(geneInfo, file = "Gene Modules Merged IP DeSeq_Minus.csv")

# Save Module Hub Genes
# hubs = chooseTopHubInEachModule(GeneCounts, mergedColors)
# write.csv(hubs, file = "Hubs Merged IP DeSeq_Minus.csv",  row.names=FALSE)

# Save Module Eigengenes Merged
# write.csv(mergedMEs, file = "MEs Merged IP DeSeq_Minus.csv",  row.names=FALSE)
