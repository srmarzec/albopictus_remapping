####
###
##
# NOTE: This is an incomplete script that so far gives KEGG IDs and gene names if given a list of gene IDs
##
###
####

# KEGG analysis

# BiocManager::install("KEGGREST")

# BiocManager::install("mygene")

library(KEGGREST)

gene_list <- read.csv("../misc/11D_DvND_test.txt", header = F)

# Make a list of all the ncbi to albopictus gene ids
convs <- keggConv("ncbi-geneid", "aalb")
convs2 <- keggConv("aalb", "uniprot")
#convs3 <- keggConv("aalb", "ncbi-proteinid")

# Convert gene ids from list to something relevant to our work. Note that all our gene ids begin with "LOC". For NCBI LOC1234 is equivalent to GeneID = 1234. The "LOC"+GeneID is when orthologs have not yet been determined and a published symbol is not available. 
gene_list$ncbi_geneid <- sub("LOC","ncbi-geneid:", gene_list$V1)


# Find the matching ncbi id in the conversion list and take the name associated (the kegg id) and assign that in the kegg_id column of the gene_list dataframe
gene_list$kegg_id = names(convs)[match(gene_list$ncbi_geneid, as.character(convs))]
gene_list$uniprot_id = names(convs2)[match(gene_list$kegg_id, as.character(convs2))]
#gene_list$protein_id = names(convs3)[match(gene_list$kegg_id, as.character(convs3))] # getting the protein id does nothing for us downstream as these are not searchable in the fasta file. I think this is likely due to protein "versions" (basically all the protein ids in the fasta end with a ".#" following so I think they are alternative proteins for the genes)

write.table(gene_list$kegg_id, 
          file="../misc/11D_DvND_keggID.txt", col.names = F, row.names = F, quote = F)
###
# You can input this KEGG list at "https://www.kegg.jp/kegg/tool/map_pathway1.html" to find enriched pathways
###


# Getting gene names using mygene package to query the gene symbols (i.e. LOC#######)
library(mygene)

query(q="LOC115271008")
query(q="109413713")
dat <- queryMany(gene_list$V1, 
                 scopes="symbol", 
                 fields=c("entrezgene", "name", "taxid")) # I don't really need taxid I had it there to double check everything was albopictus
head(dat)
