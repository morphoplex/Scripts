setwd("/Users/jonbra/Dropbox/Prosjekter/Ctenophore-transcriptomics/Mnmeiopsis/DE-analysis/")

library(DESeq2)

# See http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
# and the DESeq2 vignette

#######################################################################################################################
############################ Set up stuff needed for HTSeq counts import and DESeq dataset ############################
#######################################################################################################################

#Set the directory where the HTSeq count files are located
HTSeq_antisense = "/Users/jonbra/Dropbox/Prosjekter/Ctenophore-transcriptomics/Mnmeiopsis/DE-analysis/HTSeq-files/antisense_lncRNA/"
HTSeq_intronic = "/Users/jonbra/Dropbox/Prosjekter/Ctenophore-transcriptomics/Mnmeiopsis/DE-analysis/HTSeq-files/intronic_lncRNA/"
HTSeq_lincRNA = "/Users/jonbra/Dropbox/Prosjekter/Ctenophore-transcriptomics/Mnmeiopsis/DE-analysis/HTSeq-files/lincRNA/"
HTSeq_genes_original_annot = "/Users/jonbra/Dropbox/Prosjekter/Ctenophore-transcriptomics/Mnmeiopsis/DE-analysis/HTSeq-files/original_genes/"
HTSeq_XLOCS = "/Users/jonbra/Dropbox/Prosjekter/Ctenophore-transcriptomics/Mnmeiopsis/DE-analysis/HTSeq-files/XLOC/"
HTSeq_genes_merged_fixed = "/Users/jonbra/Dropbox/Prosjekter/Ctenophore-transcriptomics/Mnmeiopsis/DE-analysis/HTSeq-files/gene_name/"

#Catch the files with the pattern "GFF"
HTSeqFiles = grep("txt", list.files(HTSeq_genes_original_annot), value=TRUE)

#Setup sample conditions in the same order of the files in sampleFiles. Not all of these conditions are used necessarily
sample = c("HTSeq_1-1-1-1_S1_L00X.txt", "HTSeq_1-2-1-3_S2_L00X.txt", "HTSeq_1-3-2-1_S3_L00X.txt", "HTSeq_1-4-2-3_S4_L00X.txt", "HTSeq_1-4B_S7.txt", "HTSeq_1-5-3-1_S5_L00X.txt", "HTSeq_1-6-3-3_S6_L00X.txt", "HTSeq_1-7-4-1_S7_L00X.txt", "HTSeq_1-8-4-3_S8_L00X.txt", "HTSeq_10-21B_S11.txt", "HTSeq_11-12M_S14.txt", "HTSeq_12-13M_S17.txt", "HTSeq_13-14M_S20.txt", "HTSeq_15-17M_S23.txt", "HTSeq_16-3T_S2.txt", "HTSeq_17-4T_S5.txt", "HTSeq_18-7T_S9.txt", "HTSeq_19-8T_S12.txt", "HTSeq_2-5B_S10.txt", "HTSeq_20-9T_S15.txt", "HTSeq_21-10T_S18.txt", "HTSeq_22-11T_S21.txt", "HTSeq_23-12T_S24.txt", "HTSeq_24-13T_S3.txt", "HTSeq_25-17T_S6.txt", "HTSeq_3-6B_S13.txt", "HTSeq_4-7B_S16.txt", "HTSeq_5-8B_S19.txt", "HTSeq_6-10B_S22.txt", "HTSeq_7-13B_S1.txt", "HTSeq_8-15B_S4.txt", "HTSeq_9-17B_S8.txt")
condition = c("Ina_aboral", "Ina_oral", "Ina_aboral", "Ina_oral", "oral_liten", "Ina_aboral", "Ina_oral", "Ina_aboral", "Ina_oral", "oral_stor", "midt_stor", "midt_stor", "midt_stor", "midt_stor", "aboral_liten", "aboral_liten", "aboral_liten", "aboral_liten", "oral_liten", "aboral_liten", "aboral_stor", "aboral_stor", "aboral_stor", "aboral_stor", "aboral_stor", "oral_liten", "oral_liten", "oral_liten", "oral_stor", "oral_stor", "oral_stor", "oral_stor")
body_part = c("aboral", "oral", "aboral", "oral", "oral", "aboral", "oral", "aboral", "oral", "oral", "midt", "midt", "midt", "midt", "aboral", "aboral", "aboral", "aboral", "oral", "aboral", "aboral", "aboral", "aboral", "aboral", "aboral", "oral", "oral", "oral", "oral", "oral", "oral", "oral")
size = c("medium", "medium", "medium", "medium", "liten", "medium", "medium", "medium", "medium", "stor", "stor", "stor", "stor", "stor", "liten", "liten", "liten", "liten", "liten", "liten", "stor", "stor", "stor", "stor", "stor", "liten", "liten", "liten", "stor", "stor", "stor", "stor")

sampleTable = data.frame(sampleName=sample, fileName=HTSeqFiles, sample=condition, size=size, body_part=body_part)

#Creates the DESeqDataSet which stores the count data for DESeq2
dds = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq_genes_original_annot, design = ~ size)
nrow(dds)
# Ut i fra PCA-plotet ser det ut til ?? v??re en del variasjon mellom replikatene, i tillegg til variasjon mellom stor og medium+liten
# tanken med en slik design er at vi tar h??yde for forskjellen i st??rrelse n??r vi tester for DE mellom topp/bunn/midt
# The design can be changed later with design(dds) = ~ something new...

# Clean up
rm(HTSeq_intronic, HTSeq_lincRNA, HTSeq_genes_original_annot, HTSeq_XLOCS, HTSeq_genes_merged_fixed,
   HTSeqFiles,
   sample, condition, body_part, size,
   sampleTable)
#######################################################################################################################


#######################################################################################################################
################################## combine the counts for original genes and lncRNAs ##################################
#######################################################################################################################

HTSeqFiles_genes = grep("txt", list.files(HTSeq_genes_original_annot), value=TRUE)
HTSeqFiles_linc = grep("txt", list.files(HTSeq_lincRNA), value=TRUE)
HTSeqFiles_intronic = grep("txt", list.files(HTSeq_intronic), value=TRUE)
HTSeqFiles_antisense = grep("txt", list.files(HTSeq_antisense), value=TRUE)

sampleTable_genes = data.frame(sampleName=sample, fileName=HTSeqFiles_genes, sample=condition, size=size, body_part=body_part)
sampleTable_linc = data.frame(sampleName=sample, fileName=HTSeqFiles_linc, sample=condition, size=size, body_part=body_part)
sampleTable_intronic = data.frame(sampleName=sample, fileName=HTSeqFiles_intronic, sample=condition, size=size, body_part=body_part)
sampleTable_antisense = data.frame(sampleName=sample, fileName=HTSeqFiles_antisense, sample=condition, size=size, body_part=body_part)

dds_genes = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable_genes, directory=HTSeq_genes_original_annot, design = ~ size + body_part)
dds_linc = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable_linc, directory=HTSeq_lincRNA, design = ~ size + body_part)
dds_intronic = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable_intronic, directory=HTSeq_intronic, design = ~ size + body_part)
dds_antisense = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable_antisense, directory=HTSeq_antisense, design = ~ size + body_part)

nrow(dds_genes)
nrow(dds_linc)
nrow(dds_intronic)
nrow(dds_antisense)

# export the raw counts
raw_counts_genes = counts(dds_genes)
raw_counts_linc = counts(dds_linc)
raw_counts_intronic = counts(dds_intronic)
raw_counts_antisense = counts(dds_antisense)

# create lists of gene names - useful for extracting genes etc.
genes_names = row.names(raw_counts_genes)
lincRNA_names = row.names(raw_counts_linc)
intronicRNA_names = row.names(raw_counts_intronic)
antisenseRNA_names = row.names(raw_counts_antisense)

# Combine the counts into one matrix
raw_counts_combined = rbind(raw_counts_genes, raw_counts_linc, raw_counts_intronic, raw_counts_antisense)


# Create the sample table (the difference from using HTSeq import function is that we use no fileName=)
coldata = data.frame(sampleName=sample, sample=condition, size=size, body_part=body_part)
rownames(coldata) <- colnames(raw_counts_combined)

# Create a new DESeqDataSet for the combined counts
dds = DESeqDataSetFromMatrix(countData =raw_counts_combined, colData = coldata, design = ~ size + body_part)
nrow(dds)

# remove temporary files
rm(HTSeqFiles_antisense, HTSeqFiles_intronic, HTSeqFiles_linc, HTSeqFiles_genes, 
   sampleTable_antisense, sampleTable_intronic, sampleTable_linc, sampleTable_genes, 
   dds_antisense, dds_intronic, dds_linc, dds_genes,
   raw_counts_antisense, raw_counts_intronic, raw_counts_linc, raw_counts_genes)
#######################################################################################################################


#######################################################################################################################
######################################### Continue with finished DESeqDataSet #########################################
#######################################################################################################################

# Remove rows of the DESeqDataSet that have no counts, or only a single count across all samples (to speed up the analyses)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

#endre design
design(dds) = ~ body_part + size # blokker for body_part? I.e. tar bort eventuelle effekter av oral/aboral?
#design(dds) = ~ size
#design(dds) = ~ size * body_part # funker ikke pga not full rank

# Run the DE-test (If you just want to transform data for plotting/visuallization etc, don't need to rung the Wald test)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
# dds = DESeq(dds) # performs all the tree commands above

############# Transformation ################

# Blind means that the transformation should be blind to the sample information specified by the design formula
# Default is Blind = TRUE

#regularized log transformation
rld = rlog(dds)
save(dds, rld, file="dds_rld_genes_merged_fixed.RData")
load("Mnemiopsis_DESeq2.RData")
############# Visualization (sample comparisons) ############# 

#Heatmaps
library("RColorBrewer")
library("gplots")

#Select the most highly expressed genes (highest mean values across rows)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

#heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
#          Rowv = FALSE, Colv = FALSE, scale="none",
#          dendrogram="none", trace="none", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = TRUE, Colv = TRUE, scale="none",
          dendrogram="both", trace="none", margin=c(10, 6))
#heatmap.2(assay(vst)[select,], col = hmcol,
#          Rowv = FALSE, Colv = FALSE, scale="none",
#          dendrogram="none", trace="none", margin=c(10, 6))

#Heatmap of sample-to-sample distances (uses all genes in the sample)
distsRL <- dist(t(assay(rld)))
# distsRL <- dist(t(counts(dds, normalized=TRUE))) # cluster using normalized counts instead of rlog transformed
mat <- as.matrix(distsRL)
# rownames(mat) <- colnames(mat) <- with(colData(dds), paste(Condition))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13), main = "Sample-to-sample distances (rlog)")

# Using the 1000 most highly expressed genes (selected by mean expression across samples):
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:2000]
distsRL <- dist(t(assay(rld)[select,]))
# mat <- as.matrix(distsRL)
hc = hclust(distsRL)
plot(hc, main="Sample-to-sample distances, 1000 top genes (rlog)")

# Selecting genes based on variance (the varFilter is in the genefilter package is especially designed for gene data)
library(genefilter)
m.var <- varFilter(assay(rld), var.func=IQR, var.cutoff=0.6, filterByQuantile=TRUE)
distsRL <- dist(t(m.var))
hc = hclust(distsRL)
plot(hc, main="Sample-to-sample distances, variable genes (rlog)")

# PCA plot
pdf("DESeq2.PCA.genes_merged_fixed.pdf")
print(plotPCA(rld, intgroup=c("body_part", "size")))
dev.off()


#######################################################################################################################
###################################################### DE testing #####################################################
#######################################################################################################################

res = results(dds)
res
save(res, file="DESeq2_results_for_goseq_test.RData")
mcols(res, use.names=TRUE)
summary(res)
table(res$padj < 0.1)
plotMA(res, ylim=c(-4,4)) # We see that there are a lot of genes significantly DE, despite low fold change. Especially for high expression. 

# filter on log change
resLFC1 <- results(dds, lfcThreshold=1) # NB. Dette er vel det samme som log2FoldChange p?? 2?? Se F1000 artikkelen. Virker ogs?? s??nn n??r jeg ser p?? tabellen i excel...
table(resLFC1$padj < 0.1)
summary(resLFC1)
plotMA(resLFC1, ylim=c(-5,5)) # We see now that filtering on log2FC of 1/-1 (which means doubling or halving, i.e. 2^1) reduces a lot
# extract genes with padj < 0.1
resLFC1Sig = resLFC1[which(resLFC1$padj < 0.1), ]
summary(resLFC1Sig)
names_resLFC1Sig=row.names(resLFC1Sig)
write.table(names_resLFC1Sig, file = "names_significant_DE_stor_vs_liten_lcf1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# upregulated
resLFC1SigUpReg = resLFC1Sig[which(resLFC1Sig$log2FoldChange > 0), ]
summary(resLFC1SigUpReg)
table(resLFC1SigUpReg$padj < 0.1)

# Liten vs. medium
resLM = results(dds, contrast = c("size", "liten", "medium"), lfcThreshold=1)
resLM
summary(resLM)
table(resLM$padj < 0.1)

#vary the alpha:
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)

hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")


### Extract the different gene sets
res_genes = res[which(row.names(res) %in% genes_names), ]
res_linc = res[which(row.names(res) %in% lincRNA_names), ]
res_intronic = res[which(row.names(res) %in% intronicRNA_names), ]
res_antisense = res[which(row.names(res) %in% antisenseRNA_names), ]

res_antisenseSig = res_antisense[which(res_antisense$padj < 0.1), ]
summary(res_antisenseSig)
names_antisese_Sig = row.names(res_antisenseSig)
write.table(names_antisese_Sig, file = "names_significant_DE_antisenselncRNAs.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#export the DE genes to use in GOSeq
DEgenes=as.integer(p.adjust(res$pvalue[res$log2FoldChange!=0], method="BH")<.1) # bør kunne la vær å kjøre p.adjust, men heller filtrere direkte p?? padj?
names(DEgenes)=row.names(res[res$log2FoldChange!=0,])
save(DEgenes, file="DEgene_list_for_GOseq.R")

# extract results from different conditions
colData(dds)
res.size = results(dds, contrast=c("size", "liten", "stor"))
summary(res.size)
res.part = results(dds, contrast=c("body_part", "oral", "aboral")) # samme som res uten contrast
res.part = res.part[order(res.part$padj),]
res.part
summary(res.part)
res.m.o = results(dds, contrast=c("body_part", "oral", "midt"))
summary(res.m.o)
res.m.a = results(dds, contrast=c("body_part", "aboral", "midt"))
summary(res.m.a)

### Plotting counts for individual genes
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("body_part"))

library(ggplot2)
data <- plotCounts(dds, gene=topGene, intgroup=c("body_part","size"), returnData=TRUE)
ggplot(data, aes(x=body_part, y=count, color=size)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=3)

ggplot(data, aes(x=body_part, y=count, color=size, group=size)) +
  scale_y_log10() + geom_point(size=3) + geom_line()



library("genefilter")
library(pheatmap)
topVarGenes <- head(order(-rowVars(assay(rld))),20) # extract the 20 most variable genes
mat <- assay(rld)[ topVarGenes, ]

namesLFC1Sig = row.names(resLFC1Sig[resLFC1Sig$padj < 0.1,])
matLFC1Sig = assay(rld)[namesLFC1Sig,]

mat <- matLFC1Sig - rowMeans(matLFC1Sig) # To look at the amount by which each gene deviates in a specific sample from the gene???s average across all samples. Hence, we center each genes??? values across samples,
df <- as.data.frame(colData(rld)[,c("size","body_part")])
pheatmap(mat, annotation_col=df)

# Alternativt heatmap
library("RColorBrewer")
library("gplots")

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(matLFC1Sig, col = hmcol, trace = "none", Colv = F, Rowv = T, scale = "row") # clustre kun rader
# kun stor og liten samples:
heatmap.2(matLFC1Sig[,c("HTSeq_1-4B_S7.txt", "HTSeq_10-21B_S11.txt", "HTSeq_11-12M_S14.txt", "HTSeq_12-13M_S17.txt", "HTSeq_13-14M_S20.txt", "HTSeq_15-17M_S23.txt", "HTSeq_16-3T_S2.txt", "HTSeq_17-4T_S5.txt", "HTSeq_18-7T_S9.txt", "HTSeq_19-8T_S12.txt", "HTSeq_2-5B_S10.txt", "HTSeq_20-9T_S15.txt", "HTSeq_21-10T_S18.txt", "HTSeq_22-11T_S21.txt", "HTSeq_23-12T_S24.txt", "HTSeq_24-13T_S3.txt", "HTSeq_25-17T_S6.txt", "HTSeq_3-6B_S13.txt", "HTSeq_4-7B_S16.txt", "HTSeq_5-8B_S19.txt", "HTSeq_6-10B_S22.txt", "HTSeq_7-13B_S1.txt", "HTSeq_8-15B_S4.txt", "HTSeq_9-17B_S8.txt")], col = hmcol, trace = "none", Colv = T, Rowv = T, scale = "row")
heatmap.2(matLFC1Sig, col = hmcol, trace = "none", Colv = T, Rowv = T, scale = "row") # clustre b??de rader og kolonner
heatmap.2(scale(matLFC1Sig), col = hmcol, trace = "none", Colv = F, scale = "none") # scaling before clustering
heatmap.2(scale(matLFC1Sig), col = greenred, trace = "none", Colv = F, scale = "none") # color greenred

# time course 
ddsTC = DESeq(dds, test="LRT", reduced = ~ size)
resTC = results(ddsTC)
resTC$symbol = mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)

data <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup=c("size","body_part"), returnData=TRUE)
                              ggplot(data, aes(x=body_part, y=count, color=size, group=body_part)) +
                              geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10()
