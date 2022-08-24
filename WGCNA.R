

# Display the current working directory
setwd("/Users/brianmckinley/desktop/wgcna_test/")
library(WGCNA)
options(stringsAsFactors = FALSE);
rawData = read.csv("4.1_atlas_GRN_15XDE_20TPM_genes.csv");
# Take a quick look at what is in the data set:
dim(rawData);
names(rawData);

datExpr0 = as.data.frame(t(rawData[, -c(1)]));
names(datExpr0) = rawData$Gene_ID;
rownames(datExpr0) = names(rawData)[-c(1)];

powers2 = c(seq(from = 1, to = 20, by = 1))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, 
                        powerVector = powers2, 
                        RsquaredCut = 0.7, 
                        verbose = 5, 
                        moreNetworkConcepts = TRUE, 
                        networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,1));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers2,cex=cex1,col="red");
abline(h=0.75,col="red")


# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers2, cex=cex1,col="red")

softPower = 14

adjacency = adjacency(datExpr0, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency,TOMType = "signed");


cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("WGCNA_TOM", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = names(datExpr0));




