setwd("/Users/jonbra/Dropbox/Prosjekter/Ctenophore-transcriptomics/Mnmeiopsis/DE-analysis/")

library(GenomicFeatures)
library(rtracklayer)

# txdb = makeTxDbFromGFF("../merged_fixed.gtf", format = ("gtf"))

# I have removed the genes without an "old name". 
# I.e. extracted only the lines with "gene_name" and replaced that with "gene_id"
# This is because the transcriptsBy function below uses the "gene_id" feature
txdb = makeTxDbFromGFF("../merged_fixed_for_goseq.gtf", format = ("gtf")) 

# tx_import = import("../merged_fixed.gtf", format = ("gtf"))
# txdb = makeTxDbFromGRanges(tx_import)

#Create a vector named by gene ID containing the median transcript length of each gene
txsByGene=transcriptsBy(txdb,"gene")
lengthData=median(width(txsByGene))

# Read the GO-associations file. This file was created on the basis of the Blast2Go results. Can also be used for Ontologizer
assoc = read.table(file="GO_associaton_for_R.txt", sep = " ", header=F)

#Read in the DE-genes. This list was created in the DESeq2 analysis. Will differ depending on which stages are compared. 
#See the script Mnemiopsis_SESeq2.R
load("DEgene_list_for_GOseq.R")


#### GO testing
library(goseq)
# remove NAs
DEgenes_new <- na.omit(DEgenes)
# make lengthData same length as DE_genes, i.e. contain the same genes. lengthData contains all the genes in the gtf-file,
# but for instance in the DESeq2 step we removed genes with no counts. 
DElength = lengthData[names(lengthData) %in% names(DEgenes_new)]

# Weight each gene according to its length
pwf = nullp(DEgenes_new, bias.data = DElength)
head(pwf)
nrow(pwf) # the pwf contains all the genes in DElength. I guess they are supplied from the DElength object. And probably these are used as the reference set.

GO.wall=goseq(pwf, gene2cat = assoc, method='Wallenius')
head(GO.wall)

#Using random sampling to check consistency
GO.samp=goseq(pwf,gene2cat = assoc ,method="Sampling",repcnt=10000)
head(GO.samp)

# plot the p-values against each other to see differences (should not be much)
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]), 
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
abline(0,1,col=3,lty=2) # There are very little difference

# testing for length bias (length bias is accounted for in the GO.wall above)
GO.nobias=goseq(pwf,gene2cat = assoc, method="Hypergeometric")
head(GO.nobias)
plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.nobias[,1],GO.wall[,1]),2]), 
     xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)", 
     xlim=c(-3,0), ylim=c(-3,0))
abline(0,1,col=3,lty=2)
# There are some differences when not accounting for length bias, but not very much


# Getting GO categories over enriched using a .1 FDR cutoff (BH-adjusted)
enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.1]
enriched.GO

# Loading the gene ontology database
library(GO.db)
for(go in enriched.GO[1:7]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

# Getting GO categories under enriched using a .1 FDR cutoff (BH-adjusted)
U_enriched.GO = GO.wall$category[p.adjust(GO.wall$under_represented_pvalue, method="BH")<.1]
U_enriched.GO

for(go in U_enriched.GO[1:7]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}


###### Alternative method to obtain a reference/background set for enrichment test ###### 

# https://support.bioconductor.org/p/57884/
#I guess the most important thing in an enrichment analysis is to choose a proper "background" set of genes. 
# ("measured genes" in the GOSeq parlance). In order to choose an appropriate set, you could adapt the code below. 
# It first extracts the differentially expressed genes and then uses the function matchit to compile a background 
# set that has a comparable distribution (My comment: I guess they mean distribution of counts). 
# This way, you can be sure that your background does not have a dramatically different distribution from your DE genes.

# Best wishes, Bernd

### extract differentially expressed genes
load(file="DESeq2_results_for_goseq_test.RData")
genes <- rownames(subset(res, padj <0.1))

#Now, we generate a background set of genes matched for
#expression strength, avoiding potential biases.

library(MatchIt)
backM <- c()

## prepare data frame for matching, sign indicates wheather
## the gene is differentially expressed or not
df <- data.frame( sign=as.numeric( rownames(res) %in% as.character(genes)), res["baseMean"])
df <- data.frame( sign=as.numeric( rownames(res) %in% as.character(DEgenes)), res["baseMean"])

df$baseMean <- round( df$baseMean, 0)

## repeat matching multiple times since
## each differentially expressed gene is
## matched by exactly one non-expressed gene

for( i in 1:3 ){
  mm <- matchit( sign ~ baseMean, df, 
                 method="nearest", distance="mahalanobis")
  backM <- c(backM, mm$match.matrix[,1]) 
  df <- df[which(!rownames(df) %in% backM),] 
}

backM <- unique( na.exclude(backM) )

## total number of matched genes
length(backM )

## no DE genes in Background:
intersect(backM, genes)

save(backM, file="reference_set_for_goseq.RData")
#----------------------------------------------------------
###### Compare distributions of background and significant genes
#----------------------------------------------------------

#We check whether our background has the same distribution of expression strength.

pdf("matching.pdf" , height = 5 , width = 5)
library(geneplotter)
multidensity( list( all= res[,"baseMean"] , fore=res[rownames(res) %in% genes, "baseMean"], back=res[rownames(res) %in% backM, "baseMean"]), xlab="mean counts", xlim=c(0, 150))
dev.off()
