# POGO-make-phyloseq.R
# Subset of Knight Inlet data
# author: Libby Natola
# begun: 30 April 2024

# make phyloseq object of 

setwd("/mnt/Genomics/Working/libby.natola/home/libby.natola/projects/POGO/")

# load libraries ----
library(tidyverse)
library(phyloseq)
library(data.table)
library(knitr)
library(microbiome)
library(ranacapa)
library(qualpalr) # for color palettes
library(lubridate)  # For date manipulation
library(vegan) # for beta dispersion
library(plotly) # for interactive plots
library(taxize) # to clean taxa tables
library(microViz) # to clean taxa tables and for heatmaps
library(patchwork) # for heatmaps
library(ggpubr) # for heatmaps
library(ggrepel) # for jitterin plot labels so they are readable
library(ggnewscale) # for multiple color scales



# data processing ----
## load in data ----
# sequence table
seqtab.COI <- fread("evan-files/POGO/sequence_table.CO1.merged.w_ASV_names.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.COI) <- seqtab.COI$row_names
seqtab.COI$row_names <- NULL
seqtab.COI <- t(seqtab.COI)
seqtab.COI <- seqtab.COI[, !colnames(seqtab.COI) %in% "Undetermined"] # remove undetermined sample 
seqtab.COI <- as.matrix(seqtab.COI)
# dim(seqtab.COI) [1] 3054   81

# rename the negatives so they'll match the metadata
#colnames(seqtab.COI) <- gsub("COI_NEG", "NEG", colnames(seqtab.COI))

#taxonomy table
taxa.COI <- fread("evan-files/POGO/taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.LCA+best_hit.txt", sep="\t", header=T, colClasses = c("#Query"="character"), data.table=FALSE)
row.names(taxa.COI) <- taxa.COI$`#Query`
taxa.COI$`#Query` <- NULL
taxa.COI <- as.matrix(taxa.COI[,3:9]) #separate taxonomy table from LCA info
colnames(taxa.COI) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxa.COI <- as.matrix(taxa.COI)
#dim(taxa.COI) [1] 528    7

#metadata 
metadata <- fread("Knight_IPCA_eDNA_Metadata_2022 - Sheet1.csv", sep=",", header=T, colClasses = c("Sample_ID"="character"), data.table=FALSE)
colnames(metadata)[which(colnames(metadata) == "Sample_ID")] <- "SampleID"
metadata$SampleID <- gsub("_", "-", metadata$SampleID)
row.names(metadata) <- metadata$SampleID


# check that all the metadata samples are in the seqtab data
metadata$SampleID[!metadata$SampleID %in% colnames(seqtab.COI)]

# remove metadata samples that are negatives
metadata <- metadata[!metadata$SampleType == "Negative",]

# check if all the samples i want from the seqtab are in the metadata. There are a bunch of QMI samples that's ok
notinmeta <- setdiff(colnames(seqtab.COI), row.names(metadata))
seqtab.COI <- seqtab.COI[,!(colnames(seqtab.COI) %in% notinmeta)]

# remove samples that aren't in the metadata 
notinasv <- setdiff(row.names(metadata), colnames(seqtab.COI))
metadata_ps <- subset(metadata, !(SampleID %in% notinasv))



#combine into phyloseq object
project_data.COI <- phyloseq(otu_table(seqtab.COI, taxa_are_rows=T), #taxa_are_rows is a parameter that sets whether phyloseq will look for taxa in the rows or columns of the input matrix. make sure it is set correctly.
                             sample_data(metadata_ps), 
                             tax_table(taxa.COI))


sample_data(project_data.COI)$reads_sample <- readcount(project_data.COI)
print(sample_data(project_data.COI)[,c("SampleID", "reads_sample")])

#remove zero count ASVs - went from 9912 to 6437 taxa (390 samples)
tokeep <- names(which(taxa_sums(project_data.COI) > 0))
project_data.COI <- prune_taxa(tokeep, project_data.COI)

### write output to disk 
saveRDS(project_data.COI, "project_data.RDS")
#write.table(data.frame("row_names"=rownames(asv.tab.norm),asv.tab.norm),"../r_output/50WS_ASV_table_normalized.txt", row.names=FALSE, quote=F, sep="\t")
