# Script to count up the reads assigned to genes by HTSeq for each of the samples

#if (!require("readbulk")) install.packages("readbulk")
#if (!require("tidyverse")) install.packages("tidyverse")
#tidyverse can be used to reformat tables

library(readbulk)
library(tidyverse)


# Read in all the files at once
# Now all the data in the HTSeq files are in a table called "raw_data"
raw_data <- read_bulk(directory = "/Users/cottonellezhou/OneDrive - Georgetown University/Differential Expression Analysis/Data", fun = read.table, header = F)

# Spread out the data into long format so that we can have relevant column names
# Before the data from the different samples were simply stacked on each other
# spread_dat <- spread(tablename, thecolumnnheaderofthecolumncontainingthethingsyouwanttouseasthenewcolumnheaders, thevaluesyouwanttobespreadingoutforeachofthenewcolumns)
# tidyr::spread		Spread a key-value pair across multiple columns
# Hence, must do library(tidyverse) before this command:
spread_dat <- spread(raw_data, File, V2)

# Remove the first 5 rows of reads as these are not assigned to genes
# removing rows: tablename <newtablename[c(-rownumbers)]
spread_dat <- spread_dat[c(-1,-2,-3,-4,-5),]

# get the sums for the different columns (Files/samples) knowing that this will be the number of reads of assigned to a gene
# setting variable called col_dat equal to the sums for the different columns
# colSums: Form row and column sums and means for numeric arrays (or data frames).

col_dat <- colSums(spread_dat[,-1])

# Make this into a nice little dataframe instead of a named list
col_dat <- as.data.frame(col_dat)

# Write this out as a csv so we can open it in excel later and reference the number of assigned reads (note that I am putting the full directory path here since we did not set the working directory above)
setwd("~/OneDrive - Georgetown University/Differential Expression Analysis/Script")
write.csv(col_dat, file = "../misc/htseqCount_ReadsAssignedGene.csv")
