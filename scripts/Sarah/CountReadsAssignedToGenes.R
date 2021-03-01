# Script to count up the reads assigned to genes by HTSeq for each of the samples

library(readbulk)
library(tidyverse)


# Read in all the files at once
raw_data <- read_bulk(directory = "/Users/sarah/OneDrive/Documents/Mosquito/Remapping_RNASeq/htseq-count/data", fun = read.table, header = F)

# Spread out the data into long format so that we can have relevant column names
spread_dat <- spread(raw_data, File, V2)

# Remove the first 5 rows of reads as these are not assigned to genes
spread_dat <- spread_dat[c(-1,-2,-3,-4,-5),]

#get the sums for the different columns (Files/samples) knowing that this will be the number of reads of assigned to a gene
col_dat <- colSums(spread_dat[,-1])

# Make this into a nice little dataframe instead of a named list
col_dat <- as.data.frame(col_dat)

# Write this out as a csv so we can open it in excel later and reference the number of assigned reads (note that I am putting the full directory path here since we did not set the working directory above)
write.csv(col_dat, file = "/Users/sarah/OneDrive/Documents/Mosquito/Remapping_RNASeq/htseq-count/misc/htseqCount_ReadsAssignedGene.csv")
