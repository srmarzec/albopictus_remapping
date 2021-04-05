# KEGG analysis
# The pathway enrichment analysis for this I modified from "https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html"

 BiocManager::install("KEGGREST")

# Load libraries

library(KEGGREST)

# Set some variables
setwd("~/Downloads/misc")

genes_with_pvalues <- "72h_DvND_LFCshrink_padj.txt" # with at least columns "gene", "padj"

keggGeneID_output <- "72h_DvND_keggID.txt"

# Read in data
gene_list <- read.csv(genes_with_pvalues, header = T)

# Make a list of all the ncbi to albopictus gene ids -- do we not do 2 and 3?
convs <- keggConv("ncbi-geneid", "aalb")
#convs2 <- keggConv("aalb", "uniprot")
#convs3 <- keggConv("aalb", "ncbi-proteinid")

# Convert gene ids from list to something relevant to our work. Note that all our gene ids begin with "LOC". For NCBI LOC1234 is equivalent to GeneID = 1234. The "LOC"+GeneID is when orthologs have not yet been determined and a published symbol is not available. 
gene_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", gene_list$gene)


# Find the matching ncbi id in the conversion list and take the name associated (the kegg id) and assign that in the kegg_id column of the gene_list dataframe- same thing here?
gene_list$kegg_id = names(convs)[match(gene_list$ncbi_geneid, as.character(convs))]
#gene_list$uniprot_id = names(convs2)[match(gene_list$kegg_id, as.character(convs2))]
#gene_list$protein_id = names(convs3)[match(gene_list$kegg_id, as.character(convs3))] # getting the protein id does nothing for us downstream as these are not searchable in the fasta file. I think this is likely due to protein "versions" (basically all the protein ids in the fasta end with a ".#" following so I think they are alternative proteins for the genes)

# If you want to write out the KEGG ID list and do this online
write.table(gene_list$kegg_id, 
          file=keggGeneID_output, col.names = F, row.names = F, quote = F)
###
# You can input this KEGG list at "https://www.kegg.jp/kegg/tool/map_pathway1.html" to find enriched pathways (do we do this step?)
###


########
######
####
#Trying to automate this process with KEGG pathway enrichment 
all_genes_list <- read.csv(file="../misc/DESeq_results_embryo.csv")

all_genes_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", all_genes_list$geneID)
all_genes_list$kegg_id = names(convs)[match(all_genes_list$ncbi_geneid, as.character(convs))]


# Get the pathways list from KEGG
pathways.list <- keggList("pathway", "aalb")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)})
head(genes.by.pathway)

geneList <- all_genes_list$X72h_padj
names(geneList) <- sub("aalb:","", gene_list$kegg_id) # get rid of the beginning "aalb:" since the gene list we brought from kegg doesn't have this
head(geneList)


pathway_pval <- data.frame()

#Here I get the error message: Error in wilcox.test.default(scores.in.pathway, scores.not.in.pathway,  : 
  #not enough (non-missing) 'x' observations. proceeded with the script at line 94- not sure if this changes anything?
                                           
for (pathway in 1:length(genes.by.pathway)){
      pathway.genes <- genes.by.pathway[[pathway]]
      if (!is.na(pathway.genes)){
        list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
        list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
        scores.in.pathway <- geneList[list.genes.in.pathway]
        scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
        if (length(scores.in.pathway) > 0){
            p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
        } else{
            p.value <- NA
        }
        new_row <- c(names(genes.by.pathway[pathway]), p.value, length(list.genes.in.pathway))
        pathway_pval <- rbind(pathway_pval, new_row)
      }
    }

colnames(pathway_pval) <- c("pathwayCode", "pval", "annotated")
pathway_pval <- pathway_pval[complete.cases(pathway_pval),]

pathway_pval$pathwayName = pathways.list[match(pathway_pval$pathwayCode, sub("path:","", names(pathways.list)))]

head(pathway_pval)
pathway_pval$pval <- as.numeric(pathway_pval$pval)

pathway_pval <- pathway_pval[order(pathway_pval$pval),]
head(pathway_pval)

# Write out a csv with these data
write.csv(pathway_pval, 
          file="../output/larva_72h_keggPathwayEnrichment.csv", row.names = F)  
                                           
#Rerunning with the second set of data                                           
# KEGG analysis
# The pathway enrichment analysis for this I modified from "https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html"

 BiocManager::install("KEGGREST")

# Load libraries

library(KEGGREST)

# Set some variables
setwd("~/Downloads/misc")

genes_with_pvalues <- "135h_DvND_LFCshrink_padj.txt" # with at least columns "gene", "padj"

keggGeneID_output <- "135h_DvND_keggID.txt"

# Read in data
gene_list <- read.csv(genes_with_pvalues, header = T)

# Make a list of all the ncbi to albopictus gene ids -- do we not do 2 and 3?
convs <- keggConv("ncbi-geneid", "aalb")
#convs2 <- keggConv("aalb", "uniprot")
#convs3 <- keggConv("aalb", "ncbi-proteinid")

# Convert gene ids from list to something relevant to our work. Note that all our gene ids begin with "LOC". For NCBI LOC1234 is equivalent to GeneID = 1234. The "LOC"+GeneID is when orthologs have not yet been determined and a published symbol is not available. 
gene_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", gene_list$gene)


# Find the matching ncbi id in the conversion list and take the name associated (the kegg id) and assign that in the kegg_id column of the gene_list dataframe- same thing here?
gene_list$kegg_id = names(convs)[match(gene_list$ncbi_geneid, as.character(convs))]
#gene_list$uniprot_id = names(convs2)[match(gene_list$kegg_id, as.character(convs2))]
#gene_list$protein_id = names(convs3)[match(gene_list$kegg_id, as.character(convs3))] # getting the protein id does nothing for us downstream as these are not searchable in the fasta file. I think this is likely due to protein "versions" (basically all the protein ids in the fasta end with a ".#" following so I think they are alternative proteins for the genes)

# If you want to write out the KEGG ID list and do this online
write.table(gene_list$kegg_id, 
          file=keggGeneID_output, col.names = F, row.names = F, quote = F)
###
# You can input this KEGG list at "https://www.kegg.jp/kegg/tool/map_pathway1.html" to find enriched pathways (do we do this step?)
###


########
######
####
#Trying to automate this process with KEGG pathway enrichment 
all_genes_list <- read.csv(file="../misc/DESeq_results_embryo.csv")

all_genes_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", all_genes_list$geneID)
all_genes_list$kegg_id = names(convs)[match(all_genes_list$ncbi_geneid, as.character(convs))]


# Get the pathways list from KEGG
pathways.list <- keggList("pathway", "aalb")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)})
head(genes.by.pathway)

geneList <- all_genes_list$X135h_padj
names(geneList) <- sub("aalb:","", gene_list$kegg_id) # get rid of the beginning "aalb:" since the gene list we brought from kegg doesn't have this
head(geneList)


pathway_pval <- data.frame()

#Here I get the error message: Error in wilcox.test.default(scores.in.pathway, scores.not.in.pathway,  : 
  #not enough (non-missing) 'x' observations. proceeded with the script at line 94- not sure if this changes anything?
                                           
for (pathway in 1:length(genes.by.pathway)){
      pathway.genes <- genes.by.pathway[[pathway]]
      if (!is.na(pathway.genes)){
        list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
        list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
        scores.in.pathway <- geneList[list.genes.in.pathway]
        scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
        if (length(scores.in.pathway) > 0){
            p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
        } else{
            p.value <- NA
        }
        new_row <- c(names(genes.by.pathway[pathway]), p.value, length(list.genes.in.pathway))
        pathway_pval <- rbind(pathway_pval, new_row)
      }
    }

colnames(pathway_pval) <- c("pathwayCode", "pval", "annotated")
pathway_pval <- pathway_pval[complete.cases(pathway_pval),]

pathway_pval$pathwayName = pathways.list[match(pathway_pval$pathwayCode, sub("path:","", names(pathways.list)))]

head(pathway_pval)
pathway_pval$pval <- as.numeric(pathway_pval$pval)

pathway_pval <- pathway_pval[order(pathway_pval$pval),]
head(pathway_pval)

# Write out a csv with these data
write.csv(pathway_pval, 
          file="../output/larva_135h_keggPathwayEnrichment.csv", row.names = F)                                            
