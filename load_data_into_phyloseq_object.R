####code for getting your data (from dada2 pipeline) into a phyloseq object####
#author: Evan Morien
#last modified: April 28th, 2025

#IMPORTANT NOTE: for purposes of clarity, there are example files for each of these methods of loading in data located in the lab git repo. they are stored in the subfolder "git:parfreylab/lab_scripts/example_files". these files have the format of the output files from dada2 or MED/QIIME pipelines.

#### For data created using the dada2 pipeline: load either .RDS (r object) or .txt formatted sequence table, taxonomy table, and metadata####
#loading from .RDS
phyloseq_object_name <- readRDS(name_of_file_where_you_stored_phyloseq_object.RDS) #modify for whatever file you need to load. NB this works for whatever R objects you like, you can save and load any R object, preserving all formatting, class attributes, etc and load it back in on another computer this way.

####loading sequence/taxonomy tables from tab-delimited .txt ####
setwd("/path/to/example_files") #to run this code, simply set the correct path for the example files that you cloned/downloaded from our lab's github repo (see line 5 of this script)
library(phyloseq)
library(data.table)
library(tidyverse)
seqtab.nosingletons.nochim <- fread("sequence_table.example.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with the row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix

taxatable <- fread("taxonomy_table.example.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(taxatable) <- taxatable[,1] #set row names
taxatable <- taxatable[,-1] #remove column with the row names in it
taxatable <- as.matrix(taxatable) #cast the object as a matrix

#load sample data (from tab-delimited .txt)
rawmetadata <- fread("mapping_file.example.txt", sep="\t", header=T, colClasses = c("sampleID"="character"), data.table=FALSE) #change "sampleID" to the header of the column where your sample names are stored
row.names(rawmetadata) <- rawmetadata[,1] #set row names #change the index to the column with the row names in it
rawmetadata <- rawmetadata[,-1] #remove column with the row names in it

#make a note of which samples are present/absent in the metadata vs the sequence table (samples without a corresponding entry in both of these will not be included in the final phyloseq object, so correct any mislabeled samples now (the matching must be exact)
notinmeta <- setdiff(row.names(seqtab.nosingletons.nochim), row.names(rawmetadata))
notinraw <- setdiff(row.names(rawmetadata), row.names(seqtab.nosingletons.nochim))

#combine into phyloseq object
project_data <- phyloseq(otu_table(seqtab.nosingletons.nochim, taxa_are_rows=FALSE), #taxa_are_rows is a parameter that sets whether phyloseq will look for taxa in the rows or columns of the input matrix. make sure it is set correctly.
                          sample_data(rawmetadata), 
                          tax_table(taxatable))


#### saving phyloseq object ####
# at this point I recommend saving your full unfiltered dataset in phyloseq format as a .RDS, so that you can pick up the analysis from this point easily if you decide to change your downstream filtering criteria later on
saveRDS(project_data, "my_project.full_dataset.phyloseq_format.RDS")

#### OPTIONAL but GENERALLY RECOMMENDED UNLESS THERE IS A GOOD REASON NOT TO DO THIS: filtering the complete phyloseq dataset ####
# Remove samples with less than N reads (N=100 in example. adjust per experiment.)
project_data.filt <- prune_samples(sample_sums(project_data) >= 100, project_data) #i generally err on the side of caution for 18s experiments, where low read counts can reflect low diversity/low numbers of 18s organisms in the sample, while not strictly reflecting "sample failure". for 16s, a higher threshold is better. successful samples will have more than 1000 reads, but use a plot of the sorted sample sums to find what looks like an appropriate data-informed cutoff, if possible.

# OPTIONAL: Remove OTUs with less than N total reads. (N = 50 for example. adjust per experiment)
project_data.filt <- prune_taxa(taxa_sums(project_data) >= 50, project_data) #again, discretion is necessary here. in the dada2 pipeline, this is not necessary becasue a similar filter is applied already in the pipeline. for MED/QIIME, user's discretion.

# Set aside unassigned taxa #with dada2, there may not be any unassigned taxa as dada2's RDP classifier usually ends up classifying everything. You may want to adjust this command to set aside ASVs that have no assignment beyond Rank1 instead of no assignment at Rank1.
project_data.unassigned <- project_data.filt %>%
  subset_taxa(Rank1 == "Unassigned")
# Remove unassigned taxa
project_data.filt <- project_data.filt %>%
  subset_taxa(Rank1 != "Unassigned")

# zero out counts of 2 or less from OTU table
# this is a de-noising procedure applicable to both dada2 and MED/QIIME produced datasets which minimizes cross-well contamination
otu <- as.data.frame(otu_table(project_data))
otu_table(project_data)[otu <= 2] <- 0

# additional filtering steps, if not already done before loading into R
# 16S ONLY: Remove mitochondrial and chloroplast OTUs
# IMPORTANT: make sure that the filter terms will work with your taxonomy strings, and ranks. mitochondria and chloroplast may be at different ranks in whatever database you are using. check your taxonomy table manually!
project_data <- project_data %>%
  subset_taxa(Rank5 != "Mitochondria") %>%
  subset_taxa(Rank4 != "Chloroplast")

# 18S (and optional for 16S): Remove unwanted clades
# you can chain as many of these subset_taxa calls as you like into this single command using the R pipe (%>%)
project_data <- project_data %>%
  subset_taxa(RankX != "UNWANTED_CLADE") %>%
  subset_taxa(RankY != "UNWANTED_CLADE")

# modify Rank labels in taxa table (check the colnames of the tax_table(project_data) object to see if you want to change them)
colnames(tax_table(project_data)) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Accession") #SILVA ALL 16s, SILVA v132 18s
colnames(tax_table(project_data)) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Rank8", "Accession") #SILVA 128 18s, PR2 2019

#### accessing data within phyloseq objects ####
#this can be a bit tricky, but you just have to treat the component parts of the S4 object as data frames
#it is sometimes necessary to strip the class information out with "unclass()" before viewing the data tables with the View() function
#for viewing
View(as.data.frame(unclass(sample_data(project_data))))
View(as.data.frame(unclass(otu_table(project_data))))
View(as.data.frame(unclass(tax_table(project_data))))

#creating separate data frames for use with other programs
my_sample_data <- as.data.frame(unclass(sample_data(project_data))) #etc
