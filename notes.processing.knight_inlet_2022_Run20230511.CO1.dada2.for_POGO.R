#notes for processing knight_inlet_2022_Run20230511 CO1 data
#author: Evan Morien
#last modified: June 7th, 2022

####libraries####
#these are the libraries used in this pipeline. in a few cases, the order of library loading matters, so best not to modify it.
library(dada2)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)

####Environment Setup####
theme_set(theme_bw())
setwd("~/projects/CO1_runs/knight_inlet_2022_Run20230511/")
####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "~/projects/CO1_runs/knight_inlet_2022_Run20230511/raw_data/" # CHANGE ME to the directory containing the fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #change the pattern to match all your R1 files
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name

####fastq Quality Plots####
#if you have hundreds of samples in your dataset, it's easiest to just look at a few (10-50) fastq files
numplot <- 49 #CHANGE ME to the number of samples you want to include in the quality plot ( can be anywhere from 2 to length(fnFs) ) # plotting 49 quality plots produces a 7x7 plot, which i find useful and easy to read
a <- sample(fnFs, numplot) #randomly select a set of N samples
b <- which(fnFs %in% a) #identify the indices of those samples
plotfnFs <- fnFs[b] #use the indices to create two lists of corresponding forward and reverse files to plot
plotfnRs <- fnRs[b]
pdf("quality_plots.dada2.R1s.pdf", width = 16, height = 9) # define plot width and height. completely up to user.
  plotQualityProfile(plotfnFs) #this plots the quality profiles for each sample
dev.off()
pdf("quality_plots.dada2.R2s.pdf", width = 16, height = 9) # define plot width and height. completely up to user.
  plotQualityProfile(plotfnRs)
dev.off()

####primer removal####
FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"  ## CHANGE ME to your forward primer sequence #current primers here are the CO1-Leray primer set used by Hakai
REV <- "TANACYTCNGGRTGNCCRAARAAYCA"  ## CHANGE ME to your reverse primer sequence
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, trimRight = c(5,5), trimLeft = c(0,0), maxN = 0, multithread = 38, compress = TRUE, matchIDs=TRUE) #removing all reads with Ns in them, as we can't use them in the dada2 pipeline #here we also refer to the quality plots to determine if any read trimming needs to happen. trimming the first or last N bases at this stage can significantly improve read retention later on, if there was an issue with the run that resulted in very low quality bases at the 3' or 5' end of the reads

primerHits <- function(primer, fn) {
  # This funciton counts number of reads in which the primer is found in any orientation
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
index <- 10 #this is the index of the file we want to check for primers, within the lists "fn*s.filtN", it can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[index]]), #the index of the sample you'd like to use for this test is used here (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[index]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[index]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[index]]))

####OPTIONAL!!!!####
#REV <- REV.orients[["RevComp"]] #IMPORTANT!!! change orientation ONLY IF NECESSARY. see the online dada2 ITS workflow, section "Identify Primers" for details on when and why to do this

#### primer removal ####
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-j", 38,# -n 2 required to remove FWD and REV primer from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
#sanity check, should report zero for all orientations and read sets
index <- 10 #this is the index of the file we want to check for primers, within the lists "fn*s.cut", it can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[index]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[index]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[index]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[index]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001", full.names = TRUE)) #remember to change this so it matches ALL your file names!
cutRs <- sort(list.files(path.cut, pattern = "_R2_001", full.names = TRUE)) #remember to change this so it matches ALL your file names!

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length (definitely possible, good example is giardia). instead, defining a minimum sequence length is best.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(0,0), trimLeft = c(0, 0), trimRight = c(25,50), minLen = c(110,110),
                     maxN=c(0,0), maxEE=c(4,6), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=38) #matchIDs prevents reads from being merged if their read pair IDs don't match, minLen sets the minimum length for a read to be reatained and is an important tool for noise filtering. trimRight parameters can be figured out by viewing quality plots. I tend to trim once the average quality dips below 20-25 for an average sample
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
write.table(retained, "retained_reads.CO1.filterAndTrim_step.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)


####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf("error_rates.dada2.CO1.R1s.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
  plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()
pdf("error_rates.dada2.CO1.R2s.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
  plotErrors(errR, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()

####dereplication####
#here reads are dereplicated (only unique sequences and the number of times they are observed are retained)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####sample inference####
#see dada2 documentation for further details on this process where error rates are used to estimate whether an observed sequence is unique because of technical noise, or real biological variation
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0 merged reads, as i found out while processing a large 18s dataset. your guess is as good as mine as to why this error does or does not appear, but filtering out the samples that cause it is necessary for completion of the pipeline.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the filter and trim step. if you already did this but still got an error when merging, try the steps below
getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 50 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing
write.table(names(which(samples_to_keep == TRUE)), "samples_retained.txt", row.names=FALSE, quote=F, sep="\n")
write.table(setdiff(sample.names, names(which(samples_to_keep == TRUE))), "samples_removed.txt", row.names=FALSE, quote=F, sep="\n")

####merge paired reads####
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
pdf("length_histogram.merged_reads.length_var.pdf", width = 10, height = 8) # define plot width and height. completely up to user.
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot
dev.off()
#review after processing as sanity check. most ASVs should be w/in a few bp of the amplicon length, minus the length of the primers

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#generate counts of sample per ASV
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)

otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
removed <- length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

#Start ASV report
cat("dim(otus):", dim(otus),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("# ASVs Removed:", length(intersect(a,b)),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("dim(otus_filt):", dim(otus_filt),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=16, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
cat("dim(seqtab.nosingletons.nochim):", dim(seqtab.nosingletons.nochim),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low
cat("proprtion of chimeras:", sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons),file="ASV_report.txt",sep="\t",append=TRUE)

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100, track[,7]/track[,1]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras", "percent_retained_of_total")

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.CO1_merged.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.CO1.merged.txt", row.names=FALSE, quote=F, sep="\t")

#OPTIONAL, if needed: read in original sequence table with sequences as ASV labels.
#this is appropriate/necessary if you are loading in a sequence table produced on a remote machine, or in a separate instance of R
# seqtab.nosingletons.nochim <- fread("sequence_table.CO1.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
# row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
# seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with row names in it
# seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix
# mode(seqtab.nosingletons.nochim) <- "numeric"

#### save sequences for both ASV tables separately, and do taxonomy assignment with blast ####
#### replace the long ASV names (the actual sequences) with human-readable names ####
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.CO1.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "CO1_ASV_sequences.fasta") #save sequences with new names in fasta format

#IMPORTANT: sanity checks
colnames(seqtab.nosingletons.nochim) == ASV.seq #only proceed if this tests as true for all elements

#assign new ASV names
colnames(seqtab.nosingletons.nochim) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.CO1.merged.w_ASV_names.txt", row.names=FALSE, quote=F, sep="\t")
#at this point, you are ready to do BLAST-based taxonomy assignment

#assign taxonomy for dada2-processed CO1 amplicon data with blast using our Hakai CO1 barcode sequences
####taxonomy assignment with blast alone, using hakai barcode blast DB####
#### CO1 amplicon blast with NT ####
#assign taxonomy for dada2-processed CO1 amplicon data with blast using NCBI NT
mkdir blast_96_sim
# blast against genbank NT blast DB
blastn -task megablast -num_threads 46 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db /mnt/Genomics/Working/databases/NCBI_NT/2025-02-27/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.fasta  -out blast_96_sim/CO1_ASV_sequences.blast.out

#filter input fasta to only contain sequences with no hits in the above two blast runs
cat blast_96_sim/CO1_ASV_sequences*blast.out | cut -f1,1 | sort | uniq > blast_96_sim/blast_hit_ASVs.txt
grep -wv -f blast_96_sim/blast_hit_ASVs.txt sequence_ASVname_mapping.CO1.txt | cut -f1,1 | sort > blast_96_sim/no_blast_hit_ASVs.txt
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' blast_96_sim/no_blast_hit_ASVs.txt CO1_ASV_sequences.fasta > blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta

#blast this output with lower thresholds for similarity
mkdir blast_90_sim
# blast against genbank NT blast DB
blastn -task megablast -num_threads 46 -evalue 1e-5 -max_target_seqs 10 -perc_identity 90 -qcov_hsp_perc 50 -db /mnt/Genomics/Working/databases/NCBI_NT/2025-02-27/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta  -out blast_90_sim/CO1_ASV_sequences.blast.out

#IMPORTANT: next steps are done in R, for simplicity's sake. a custom script here would also work. this is simpler.
#we need to add taxonIDs for the customDB (adding them directly to the blast DB has not worked in the past, they don't get returned in the blast output). Using the blast output and a map of accessions to taxonIDs, we add the taxonIDs for each blast result.
library(tidyverse)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)
#combine blast results in a way that the add taxonomy and LCA scripts can handle
blastout_NCBINT <- read.delim("blast_96_sim/CO1_ASV_sequences.blast.out", sep="\t", header=F)
blastout_NCBINT_2 <- read.delim("blast_90_sim/CO1_ASV_sequences.blast.out", sep="\t", header=F)
blastout_NCBINT <- rbind(blastout_NCBINT, blastout_NCBINT_2) #join iterations for NT blast
tmp <- blastout_NCBINT[order(-blastout_NCBINT$V5),] #order descending order for percent identity
tmp <- tmp[order(tmp$V1),] #order by ASV name
blastout_NCBINT <- tmp 
#write to file
write.table(blastout_NCBINT, "blast_96_sim/CO1_ASV_sequences.combined_all.blast.out", row.names=F, col.names=F, quote=F, sep="\t")
#now quit R and continue with the remaining code in your bash shell

#avoid " keyerror: 'NA' " with python script by filtering on "NA" as a wholeword string
#explanation: occasionally, an output line from blast has an NA, which causes an error with the Simple-LCA script below. quick fix is to remove these lines (they're quite rare anyway)
grep -v -w "NA" blast_96_sim/CO1_ASV_sequences.combined_all.blast.out > blast_96_sim/tmp

#execute first step for the LCA program (adding taxonomy strings based on taxonIDs in blast output)
python2 /mnt/Genomics/Working/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim/tmp -t /mnt/Genomics/Working/databases/NCBI_taxonomy/2025-02-27/rankedlineage.dmp -m /mnt/Genomics/Working/databases/NCBI_taxonomy/2025-02-27/merged.dmp -o taxonomy
cat <(head -n 1 /mnt/Genomics/Working/programs/galaxy-tool-lca/example/example.tabular) taxonomy_tmp > tmp #here we need to add the header from the example table in order for the lca_species script to work (see below, after taxonomy string corrections)

#in this section, a series of taxonomy string modifications are made. This makes the final output more readable/easier to understand, but is optional.
#IMPORTANT: please note that the filtering criteria in the final step depend on some of this filtering (e.g. blast hits with the word "phylum" will be removed. see -fh parameter in final step) 
# because i'm replacing "unknown phylum" in these sed commands, these sequences are retained.
# if you choose not to do the replacement, the blast hits with "unknown plylum" will not be used in LCA determination unless you also change the filtering criteria set in the -fh parameter during the final step.
# also note, "unknown phylum" is present in any taxonomy where the clade's phylum is uncertain in the NCBI taxonomy system, it doesn't indicate any other kind of uncertainty about the provenance of the sequence.

#label fix for clades missing "kingdom" label
sed -i 's/unknown kingdom \/ Bacillariophyta/Eukaryota \/ Bacillariophyta/g' tmp #Bacillariophyta
sed -i 's/unknown kingdom \/ Ciliophora/Eukaryota \/ Ciliophora/g' tmp #Ciliophora
sed -i 's/unknown kingdom \/ Discosea/Eukaryota \/ Discosea/g' tmp #Discosea
sed -i 's/unknown kingdom \/ Evosea/Eukaryota \/ Evosea/g' tmp #Evosea
sed -i 's/unknown kingdom \/ Haptista/Eukaryota \/ Haptista/g' tmp #Haptista
sed -i 's/unknown kingdom \/ Rhodophyta/Eukaryota \/ Rhodophyta/g' tmp #Rhodophyta

#and for those missing kingdom + phylum labels
sed -i 's/Eukaryota \/ unknown phylum \/ Chrysophyceae/Eukaryota \/ Chrysophyceae \/ Chrysophyceae/g' tmp #Chrysophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Cryptophyceae/Eukaryota \/ Cryptophyceae \/ Cryptophyceae/g' tmp #Cryptophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Oomycota/Eukaryota \/ Oomycota \/ Oomycota/g' tmp #Oomycota
sed -i 's/Eukaryota \/ unknown phylum \/ Phaeophyceae/Eukaryota \/ Phaeophyceae \/ Phaeophyceae/g' tmp #Phaeophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Phaeophyceae/Eukaryota \/ Phaeophyceae \/ Phaeophyceae/g' tmp #Phaeophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Bigyra/Eukaryota \/ Bigyra \/ Bigyra/g' tmp #Bigyra
sed -i 's/Eukaryota \/ unknown phylum \/ Dictyochophyceae/Eukaryota \/ Dictyochophyceae \/ Dictyochophyceae/g' tmp #Dictyochophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Dinophyceae/Eukaryota \/ Dinophyceae \/ Dinophyceae/g' tmp #Dinophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Pelagophyceae/Eukaryota \/ Pelagophyceae \/ Pelagophyceae/g' tmp #Pelagophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Raphidophyceae/Eukaryota \/ Raphidophyceae \/ Raphidophyceae/g' tmp #Raphidophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Synurophyceae/Eukaryota \/ Synurophyceae \/ Synurophyceae/g' tmp #Synurophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Bolidophyceae/Eukaryota \/ Bolidophyceae \/ Bolidophyceae/g' tmp #Bolidophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Polycystinea/Eukaryota \/ Polycystinea \/ Polycystinea/g' tmp #Polycystinea
sed -i 's/Eukaryota \/ unknown phylum \/ Choanoflagellata/Eukaryota \/ Choanoflagellata \/ Choanoflagellata/g' tmp #Choanoflagellata
sed -i 's/Eukaryota \/ unknown phylum \/ Filasterea/Eukaryota \/ Filasterea \/ Filasterea/g' tmp #Filasterea

#label for those missing kindom + phylum + class labels
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Telonemida/Eukaryota \/ Telonemida \/ Telonemida \/ Telonemida/g' tmp #Telonemida
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Jakobida/Eukaryota \/ Jakobida \/ Jakobida \/ Jakobida/g' tmp #Jakobida
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Pirsoniales/Eukaryota \/ Pirsoniales \/ Pirsoniales \/ Pirsoniales/g' tmp #Pirsoniales

#execute final step for the LCA program (forming consensus taxonomies for each ASV)
mv -f tmp blast_96_sim/CO1_ASV_sequences.combined_all.blast.out #put tmp file where it belongs, add label (this is an overwrite, careful!!)
python2 /mnt/Genomics/Working/programs/galaxy-tool-lca/lca.species.py -i blast_96_sim/CO1_ASV_sequences.combined_all.blast.out -o taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.LCA+best_hit.txt -b 100 -id 90 -cov 50 -t best_hit -tid 98 -tcov 50 -fh environmental,unidentified,kingdom -flh unclassified
python2 /mnt/Genomics/Working/programs/galaxy-tool-lca/lca.species.py -i blast_96_sim/CO1_ASV_sequences.combined_all.blast.out -o taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.LCA_only.txt -b 100 -id 90 -cov 50 -t only_lca -fh environmental,unidentified,kingdom -flh unclassified

#cleanup
rm taxonomy_tmp blast_96_sim/tmp #remove redundant files