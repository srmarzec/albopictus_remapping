# This script runs DESeq (specifically starting with htseq count files)

##Install DESeq2
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install('EnhancedVolcano')
# BiocManager::install("apeglm")
# BiocManager::install("GOplot")

#Load Libraries
library(DESeq2); library(ggplot2); library(tidyverse); library(EnhancedVolcano); library(GOplot); library(pheatmap)

# Set working directory to source file location

#Choose directory with htseq-count data
directory<-("../data")

#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- grep("11d", list.files(directory), value=T) #this is pulling all files from the directory with the age specified
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- sub("_.*","", sampleFiles) #to get conditions I pull everything from before the first underscore (which I know is the diapause status)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)

# Make the DESeq dataset
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds

#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically)
dds$condition <- relevel(dds$condition, ref = "ND")

# A bit of quality control
# Look at the distribution of the count data across the samples
librarySizes <- colSums(counts(dds))

barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,1.6e+07),
        main="Barplot of library sizes")

logcounts <- log2(counts(dds) + 1)
head(logcounts)

#Is there any difference between per gene counts for each of the sample groups?
statusCol <- as.numeric(factor(dds$condition)) + 1  # make a colour vector

boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

#Adding median log counts
abline(h=median(as.matrix(logcounts)), col="blue")

#Looking at PCA of the data - Do treatments cluster together?
rlogcounts <- rlog(counts(dds))#transforming data to make it approximately homoskedastic, n < 30 so rlog is better

data_for_PCA <- t(rlogcounts)
dim(data_for_PCA)

#run PCA
pcDat <- prcomp(data_for_PCA, center = T)
pcDat_df <- data.frame('Condition' = dds$condition, pcDat$x[,1:2])

# basic plot
ggplot(pcDat_df, aes(x = PC1, y = PC2, col = Condition)) +
  geom_point() + 
  theme_bw()

# Transform normalized counts using the rlog function ("RNASEQ20_Day3_HandsOn.pdf")
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="condition")

# Hierarchial Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld) #retrieve matrix from the rld object
# compute pairwise correlation values for samples
rld_cor <- cor(rld_mat)
# plot the correlation c=values as a heatmap
pheatmap(rld_cor)

# # Look at the difference in total reads from the different samples
# dds <- estimateSizeFactors(dds)
# 
# # Estimate the dispersion (or variation) or the data
# dds <- estimateDispersions(dds)
# 
# plotDispEsts(dds)

# Differential Expression Analysis
#Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(dds)
res_LFC <- lfcShrink(dds, coef="condition_D_vs_ND", type="apeglm")
res_LFC

# Order the table by smallest p value
resOrdered <- res[order(res$pvalue),]
summary(res)


#Basic MA-plot
plotMA(res, ylim=c(-3,3))
plotMA(res_LFC, ylim=c(-3,3)) # See lots of shrinkage towards the beginning

#If you want to interactivly find genes associated with certain points
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
EnhancedVolcano(res_LFC,
                lab = rownames(res_LFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-1.5, 1.5),
                ylim = c(0,30),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)

# Can look at which of the result are significant and have a high enough log2 fold change
sig_res <- res_LFC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 0.58)

# Write out a table of these significant differentially expressed genes
write.csv(select(sig_res, gene, log2FoldChange, padj), 
            file="../misc/11D_DvND_LFCshrink_padj.txt", col.names = F, row.names = F)

# Write out just the gene names for later analysis in KEGG
write.table(sig_res %>% select(gene), 
            file="../misc/11D_DvND_test.txt", col.names = F, row.names = F, quote = F)


####################
# Running through another dataset (21days)
#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- grep("21d", list.files(directory), value=T) #this is pulling all files from the directory with the age specified
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- sub("_.*","", sampleFiles) #to get conditions I pull everything from before the first underscore (which I know is the diapause status)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)

# Make the DESeq dataset
dds_21d <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds_21d

#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds_21d)) >= 10
dds_21d <- dds_21d[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically)
dds_21d$condition <- relevel(dds_21d$condition, ref = "ND")


#Run DESeq
dds_21d <- DESeq(dds_21d)
res_21d <- results(dds_21d)
res_21d

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(dds_21d)
res_LFC_21d <- lfcShrink(dds_21d, coef="condition_D_vs_ND", type="apeglm")
res_LFC_21d


# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
EnhancedVolcano(res_LFC_21d,
                lab = rownames(res_LFC_21d),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-1.5, 1.5),
                ylim = c(0,30),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)

# Can look at which of the result are significant and have a high enough log2 fold change
sig_res_21d <- res_LFC_21d %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 0.58)

# Write out a table of these significant differentially expressed genes
write.csv(select(sig_res_21d, gene, log2FoldChange, padj), 
          file="../misc/21D_DvND_LFCshrink_padj.txt", row.names = F)


####################
# Running through another dataset (40days)
#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- grep("40d", list.files(directory), value=T) #this is pulling all files from the directory with the age specified
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- sub("_.*","", sampleFiles) #to get conditions I pull everything from before the first underscore (which I know is the diapause status)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)

# Make the DESeq dataset
dds_40d <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = directory,
                                      design = ~ condition)
dds_40d

#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds_40d)) >= 10
dds_40d <- dds_40d[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically)
dds_40d$condition <- relevel(dds_40d$condition, ref = "ND")


#Run DESeq
dds_40d <- DESeq(dds_40d)
res_40d <- results(dds_40d)
res_40d

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(dds_40d)
res_LFC_40d <- lfcShrink(dds_40d, coef="condition_D_vs_ND", type="apeglm")
res_LFC_40d


# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
EnhancedVolcano(res_LFC_40d,
                lab = rownames(res_LFC_40d),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-1.5, 1.5),
                ylim = c(0,30),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)

# Can look at which of the result are significant and have a high enough log2 fold change
sig_res_40d <- res_LFC_40d %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 0.58)

# Write out a table of these significant differentially expressed genes
write.csv(select(sig_res_40d, gene, log2FoldChange, padj), 
          file="../misc/40d_DvND_LFCshrink_padj.txt", row.names = F)



##################
# Making a venn diagram of the 3 different datasets DEGs
dat_11 <- read.csv("../misc/11d_DvND_LFCshrink_padj.txt")
dat_21 <- read.csv("../misc/21d_DvND_LFCshrink_padj.txt")
dat_40 <- read.csv("../misc/40d_DvND_LFCshrink_padj.txt")


l1 <- dat_11[, c(1,2)]
l2 <- dat_21[, c(1,2)]
l3 <- dat_40[, c(1,2)]

VennDiag <- GOVenn(l1,l2,l3, label=c('11d','21d','40d'), plot = F)
print(VennDiag$plot)
# The colors inside the labels: Red (up), Blue (down), and Yellow (contra-regulated). These colors can be changed but the order of the regulation (up, down, contra) starts at the top and goes clockwise.
