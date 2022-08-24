# Display the current working directory
setwd("/Users/brianmckinley/desktop/wgcna_test/")
library(WGCNA)
options(stringsAsFactors = FALSE);
rawData = read.csv("WGCNA_genes.csv");
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

softPower = 10

adjacency = adjacency(datExpr0, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency,TOMType = "signed");


cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("WGCNA_TOM", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = names(datExpr0));




