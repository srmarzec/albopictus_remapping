# GO-term enrichment using TopGO for 72h
# This compares a subset of genes (genes of interest; in this case, significant DEGs) to a full set of genes
# Specifically for this script, a GO-term annotations list must be provided. This list is tab-delimited and has two columns. The first column is the gene names. The second column contains the associated GO-terms. If there are multiple GO-terms, then each is list with a comma and space in between. [I made this file separately for *Aedes albopictus* and the process can be found here: "https://github.com/srmarzec/albopictus_remapping/blob/main/misc/GO_annotation.md"]
# The following script is based on the vignette ("http://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf") and these blog post ("https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html" & "http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html") 

# Install packages if necessary
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")
install.packages("formattable")

# Load libraries
library(topGO); library(Rgraphviz); library(formattable)

gcloud compute scp bananas-controller:/home/mlp134/AalbF2_GOannotations.txt /Users/mackenzieparsons/Downloads/misc
# Set up files to read in 
annotations_file <- "../misc/AalbF2_GOannotations.txt" #Contains gene IDs and GO annotations. See specifics of GO-term annotation file above
interesting_genes <- "../misc/72h_DvND_test.txt" #This file is a single column of gene IDs with no header
project_description <- "larva 72h" #Name given for creating GO data object (something like "My project" or "larva 11d")

# Read in GO annotations for genes
geneID2GO <- readMappings(file = annotations_file) 

# GO hierarchy
# topGO uses a GO hierarchy that comes with the GO.db package, so we don't need to worry about his

# Defining the genes of interest and what we will compare to ("gene universe")
# The gene universe can be all the genes that have GO annotations (for us, this is the full set Sarah created)
geneUniverse <- names(geneID2GO) 
# Then we read in the list of interesting genes (for us, this is a set of significant DEGs)
genesOfInterest <- read.table(interesting_genes,header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
# We need to say where the interesting genes appear in the gene universe set. This means creating an ojbect that lists out which genes in the gene universe are your genes of interest for topGO.
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# Create an R object that has all this data together
myGOdata <- new("topGOdata", description=project_description, ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
# The 'ontology' argument can be 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component).
# The 'description' argument has a short description of your project.
# The 'allGenes' argument specifies all the genes in the 'gene universe', and which of those are  your genes of interest.
# The 'annot' argument tells topGO how to map genes to GO annotations. 'annot'=annFUN.gene2GO means that the user provides gene-to-GO annotations, and we specify here that they are in object 'geneID2GO' (see above for how this was created).
# An optional argument is 'nodeSize': this is used to prune the GO hierarchy, eg. nodesize=10 prunes the GO hierarchy, to remove terms which have less than 10 annotated genes.

# Running the dataset name gives you a summary 
myGOdata

#graph(myGOdata)

# The list of genes of interest can be accessed using the method sigGenes():
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
#1730

#####

# Performing enrichment tests of the GO data

resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
# Selecting 'algorithm=classic' means that the GO hierarchy isn't taken into account, so each GO term is tested independently. (Usually in fact we want to take the GO hierarchy into account, so would use algorithm='weight01' for example.)
#resultKS <- runTest(myGOdata, algorithm = "classic", statistic = "ks")
#resultKS.elim <- runTest(myGOdata, algorithm = "elim", statistic = "ks")

# Typing out resultFisher will give you a summary 
resultFisher

# We can list the top ten significant results found:
allRes <- GenTable(myGOdata, 
                   weightFisher = resultFisher, 
                   orderBy = "weightFisher", ranksOf = "weightFisher", topNodes = 10)

formattable(allRes)



# We can visualise the position of the statistically significant GO terms in the GO hierarchy by using the following functions:
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
# Output description from the vingette: "The subgraph induced by the top 5 GO terms identified by the elim algorithm for scoring GO terms for enrichment. Rectangles indicate the 5 most significant terms. Rectangle color represents the relative significance, ranging from dark red (most significant) to bright yellow (least significant). For each node, some basic information is displayed. The first two lines show the GO identifier and a trimmed GO name. In the third line the raw p-value is shown. The forth line is showing the number of significant genes and the total number of genes annotated to the respective GO term."

# The following prints out the above but as a PDF. Change the prefix in order to change the folder it prints to
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "../misc/tGO", useInfo = "all", pdfSW = TRUE)

# Find gene annotated with significant GO terms
# This print all the genes from the gene universe that have the GO term
myterms = c("GO:0043044") 
mygenes <- genesInTerm(myGOdata, myterms) # get warning: Warning message:In .genesInNode(graph(object), whichGO) :Nodes not present in the graph:GO:00430
 # Output from next lines 89 to 94 says: "Term GO:0043044 genes: "
  
for (i in 1:length(myterms)) {
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- paste(mygenesforterm, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm))
}

## Repeating for 135h set 
# GO-term enrichment using TopGO for 72h
# This compares a subset of genes (genes of interest; in this case, significant DEGs) to a full set of genes
# Specifically for this script, a GO-term annotations list must be provided. This list is tab-delimited and has two columns. The first column is the gene names. The second column contains the associated GO-terms. If there are multiple GO-terms, then each is list with a comma and space in between. [I made this file separately for *Aedes albopictus* and the process can be found here: "https://github.com/srmarzec/albopictus_remapping/blob/main/misc/GO_annotation.md"]
# The following script is based on the vignette ("http://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf") and these blog post ("https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html" & "http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html") 

# Install packages if necessary
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")
install.packages("formattable")

# Load libraries
library(topGO); library(Rgraphviz); library(formattable)

gcloud compute scp bananas-controller:/home/mlp134/AalbF2_GOannotations.txt /Users/mackenzieparsons/Downloads/misc
# Set up files to read in 
annotations_file <- "../misc/AalbF2_GOannotations.txt" #Contains gene IDs and GO annotations. See specifics of GO-term annotation file above
interesting_genes <- "../misc/135h_DvND_test.txt" #This file is a single column of gene IDs with no header
project_description <- "larva 135h" #Name given for creating GO data object (something like "My project" or "larva 11d")

# Read in GO annotations for genes
geneID2GO <- readMappings(file = annotations_file) 

# GO hierarchy
# topGO uses a GO hierarchy that comes with the GO.db package, so we don't need to worry about his

# Defining the genes of interest and what we will compare to ("gene universe")
# The gene universe can be all the genes that have GO annotations (for us, this is the full set Sarah created)
geneUniverse <- names(geneID2GO) 
# Then we read in the list of interesting genes (for us, this is a set of significant DEGs)
genesOfInterest <- read.table(interesting_genes,header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
# We need to say where the interesting genes appear in the gene universe set. This means creating an ojbect that lists out which genes in the gene universe are your genes of interest for topGO.
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# Create an R object that has all this data together
myGOdata <- new("topGOdata", description=project_description, ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
# The 'ontology' argument can be 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component).
# The 'description' argument has a short description of your project.
# The 'allGenes' argument specifies all the genes in the 'gene universe', and which of those are  your genes of interest.
# The 'annot' argument tells topGO how to map genes to GO annotations. 'annot'=annFUN.gene2GO means that the user provides gene-to-GO annotations, and we specify here that they are in object 'geneID2GO' (see above for how this was created).
# An optional argument is 'nodeSize': this is used to prune the GO hierarchy, eg. nodesize=10 prunes the GO hierarchy, to remove terms which have less than 10 annotated genes.

# Running the dataset name gives you a summary 
myGOdata

#graph(myGOdata)

# The list of genes of interest can be accessed using the method sigGenes():
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
#3427

#####

# Performing enrichment tests of the GO data

resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
# Selecting 'algorithm=classic' means that the GO hierarchy isn't taken into account, so each GO term is tested independently. (Usually in fact we want to take the GO hierarchy into account, so would use algorithm='weight01' for example.)
#resultKS <- runTest(myGOdata, algorithm = "classic", statistic = "ks")
#resultKS.elim <- runTest(myGOdata, algorithm = "elim", statistic = "ks")

# Typing out resultFisher will give you a summary 
resultFisher

# We can list the top ten significant results found:
allRes <- GenTable(myGOdata, 
                   weightFisher = resultFisher, 
                   orderBy = "weightFisher", ranksOf = "weightFisher", topNodes = 10)

formattable(allRes)



# We can visualise the position of the statistically significant GO terms in the GO hierarchy by using the following functions:
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
# Output description from the vingette: "The subgraph induced by the top 5 GO terms identified by the elim algorithm for scoring GO terms for enrichment. Rectangles indicate the 5 most significant terms. Rectangle color represents the relative significance, ranging from dark red (most significant) to bright yellow (least significant). For each node, some basic information is displayed. The first two lines show the GO identifier and a trimmed GO name. In the third line the raw p-value is shown. The forth line is showing the number of significant genes and the total number of genes annotated to the respective GO term."

# The following prints out the above but as a PDF. Change the prefix in order to change the folder it prints to
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "../misc/tGO135h", useInfo = "all", pdfSW = TRUE)

# Find gene annotated with significant GO terms
# This print all the genes from the gene universe that have the GO term
myterms = c("GO:0043044") 
mygenes <- genesInTerm(myGOdata, myterms) # get warning: Warning message:In .genesInNode(graph(object), whichGO) :Nodes not present in the graph:GO:00430
 # Output from next lines 89 to 94 says: "Term GO:0043044 genes: "
  
for (i in 1:length(myterms)) {
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- paste(mygenesforterm, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm))
}
