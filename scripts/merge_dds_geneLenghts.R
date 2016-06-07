# This is for adding gene lengths to a DESeqDataSet in order to calculate fpkm values.
# You need a DESeqDataSet (called dds from now) - see the script: Mnemiopsis_DESeq2.R. I am using XLOCS here.
# You need a lengthData object. See also the goseq.R script for using original gene names.

library(GenomicFeatures)
library(rtracklayer)

txdb = makeTxDbFromGFF("merged_fixed.gtf", format = ("gtf")) 
txsByGene=transcriptsBy(txdb,"gene")
lengthData=median(width(txsByGene)) #median width of the genes

# Merging lengthData with dds:
lengthData = as.data.frame(lengthData) #convert to data frame
# Checking that the order of rows is the same (should be TRUE):
all.equal(rownames(dds), rownames(lengthData))
# Add the matrix of gene lengths to dds:
mcols(dds) <- cbind(mcols(dds), lengthData)
# Rename column lengthData to basepairs (this could probably be done in the above command):
names(mcols(dds))[names(mcols(dds)) == "lengthData"] = "basepairs"
FPKM = fpkm(dds, robust = TRUE)
