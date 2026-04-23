# a script to gather various read mapping statistics for the read-based analysis

# load libraries
library(data.table)
library(openxlsx)

# project root directory
project_root="" # e.g. "/scratch/project_number/puukkosuo"

#### gather read-use statistics
### metagenomics

## number of read pairs in samples
setwd(paste(project_root,"/metagenomics/fastqc/multiqc_data", sep = ""))

# first get the total sequences for each sample
gen_stats <- read.csv(file = "multiqc_general_stats.txt", sep = "\t", header = T)

# remove R2 samples, they carry the same information as R1 samples and we don't need duplicate information for each sample
# same number of total sequences for paired read samples
gen_stats <- gen_stats[-grep("_R2", gen_stats$Sample),]

# get the total sequences - this is the number of read pairs in the samples
tot_seqs <- as.numeric(gen_stats[,grep("total_sequences", colnames(gen_stats))])
names(tot_seqs) <- gsub("_R1", "", gen_stats$Sample)

# parse names
names(tot_seqs)  <- paste("P",unlist(lapply(names(tot_seqs), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## number of read pairs in samples in cutadapat trimmed data
setwd(paste(project_root,"/metagenomics/fastqc/cutadapt_trimmed/multiqc_data", sep = ""))

# first get the total sequences for each sample
gen_stats_trimmed <- read.csv(file = "multiqc_general_stats.txt", sep = "\t", header = T)

# remove R2 samples, they carry the same information as R1 samples and we don't need duplicate information for each sample
# same number of total sequences for paired read samples
gen_stats_trimmed <- gen_stats_trimmed[-grep("_R2", gen_stats_trimmed$Sample),]

# get the total sequences - this is the number of read pairs in the samples
tot_seqs_trimmed <- as.numeric(gen_stats_trimmed[,grep("total_sequences", colnames(gen_stats_trimmed))])
names(tot_seqs_trimmed) <- gsub("_R1", "", gen_stats_trimmed$Sample)

# parse names
names(tot_seqs_trimmed)  <- paste("P",unlist(lapply(names(tot_seqs_trimmed), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## kegg read-based
setwd(paste(project_root,"/metagenomics/kegg_diamond/combined_hits", sep = ""))

all_sample_info_files <- list.files()
total_kegg_counts <- numeric(length(all_sample_info_files))
total_ko_counts <- numeric(length(all_sample_info_files))

for(f in 1:length(all_sample_info_files)){
  # read in the mate-pair combined kegg hits 
  kegg_combined_hits <- fread(file   = all_sample_info_files[f], sep    = "\t", header = TRUE, quote  = "", fill   = FALSE, data.table = FALSE)
  
  # the number of rows is the total number of read-pairs (either R1, R2 or both) with a valid kegg gene hit
  total_kegg_counts[f] <- nrow(kegg_combined_hits)
  
  # number of read-pairs with a KO annotation, used for the KO analysis
  total_ko_counts[f] <- length(which(kegg_combined_hits$ko!=""))
}
names(total_kegg_counts) <- names(total_ko_counts) <- gsub(".txt", "", all_sample_info_files)

# parse names
names(total_kegg_counts)  <- paste("P",unlist(lapply(names(total_kegg_counts), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
names(total_ko_counts)  <- paste("P",unlist(lapply(names(total_ko_counts), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## metmarkdb read-based
setwd(paste(project_root,"/metagenomics/metmarkdb_diamond/combined_hits", sep = ""))

all_sample_info_files <- list.files()
total_metmarkdb_counts <- numeric(length(all_sample_info_files))

for(f in 1:length(all_sample_info_files)){
  # read in the mate-pair combined metmarkdb hits 
  metmarkdb_combined_hits <- fread(file   = all_sample_info_files[f], sep    = "\t", header = FALSE, quote  = "", fill   = FALSE, data.table = FALSE)
  
  # the number of rows is the total number of read-pairs (either R1, R2 or both) with a valid kegg gene hit
  total_metmarkdb_counts[f] <- nrow(metmarkdb_combined_hits)
  
  # all of these hits have a marker gene annotation
}
names(total_metmarkdb_counts) <- gsub(".txt", "", all_sample_info_files)

# parse names
names(total_metmarkdb_counts)  <- paste("P",unlist(lapply(names(total_metmarkdb_counts), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## phyloflash, read-based valid 16s counts in phyloflash output
# use the aggregated phyloflash output counts for this
setwd(paste(project_root,"/metagenomics/phyloflash/downstream", sep = ""))
load("Otu_Tax_Tables_Parsed.RData")

phyloflash_total_counts <- colSums(otu_table_all, na.rm = TRUE)

# organize all stats similarly
tot_seqs <- tot_seqs[names(phyloflash_total_counts)]
tot_seqs_trimmed <- tot_seqs_trimmed[names(phyloflash_total_counts)]
total_kegg_counts <- total_kegg_counts[names(phyloflash_total_counts)]
total_ko_counts <- total_ko_counts[names(phyloflash_total_counts)]
total_metmarkdb_counts <- total_metmarkdb_counts[names(phyloflash_total_counts)]

# put into a matrix
mg_read_use_stats <- rbind(tot_seqs, tot_seqs_trimmed, total_kegg_counts, total_ko_counts, total_metmarkdb_counts, phyloflash_total_counts)
colnames(mg_read_use_stats) <- names(phyloflash_total_counts)
rownames(mg_read_use_stats) <- c("Total_read_pairs", "Total_read_pairs_trimmed", "Total_read_kegg_hits", "Total_read_kegg_hits_with_KO", "Total_read_metmarkdb_hits", "Total_read_phyloflash_counts")
mg_read_use_stats <- data.frame(mg_read_use_stats, stringsAsFactors = F, check.names = F)

### metatranscriptomics
## number of read pairs in samples
setwd(paste(project_root,"/metatranscriptomics/fastqc/multiqc_data", sep = ""))

# first get the total sequences for each sample
gen_stats <- read.csv(file = "multiqc_general_stats.txt", sep = "\t", header = T)

# remove R2 samples, they carry the same information as R1 samples and we don't need duplicate information for each sample
# same number of total sequences for paired read samples
gen_stats <- gen_stats[-grep("_R2", gen_stats$Sample),]

# get the total sequences - this is the number of read pairs in the samples
tot_seqs <- as.numeric(gen_stats[,grep("total_sequences", colnames(gen_stats))])
names(tot_seqs) <- gsub("_R1", "", gen_stats$Sample)

# parse names
names(tot_seqs)  <- paste("P",unlist(lapply(names(tot_seqs), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## number of read pairs in samples in cutadapat trimmed data
setwd(paste(project_root,"/metatranscriptomics/fastqc/cutadapt_trimmed/multiqc_data", sep = ""))

# first get the total sequences for each sample
gen_stats_trimmed <- read.csv(file = "multiqc_general_stats.txt", sep = "\t", header = T)

# remove R2 samples, they carry the same information as R1 samples and we don't need duplicate information for each sample
# same number of total sequences for paired read samples
gen_stats_trimmed <- gen_stats_trimmed[-grep("_R2", gen_stats_trimmed$Sample),]

# get the total sequences - this is the number of read pairs in the samples
tot_seqs_trimmed <- as.numeric(gen_stats_trimmed[,grep("total_sequences", colnames(gen_stats_trimmed))])
names(tot_seqs_trimmed) <- gsub("_R1", "", gen_stats_trimmed$Sample)

# parse names
names(tot_seqs_trimmed)  <- paste("P",unlist(lapply(names(tot_seqs_trimmed), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## number of read pairs in samples in cutadapat trimmed rRNA filtered data
setwd(paste(project_root,"/metatranscriptomics/fastqc/cutadapt_trimmed/rrna_filtered/multiqc_data", sep = ""))

# first get the total sequences for each sample
gen_stats_rrna_filt <- read.csv(file = "multiqc_general_stats.txt", sep = "\t", header = T)

# remove R2 samples, they carry the same information as R1 samples and we don't need duplicate information for each sample
# same number of total sequences for paired read samples
gen_stats_rrna_filt <- gen_stats_rrna_filt[-grep("_rev", gen_stats_rrna_filt$Sample),]

# get the total sequences - this is the number of read pairs in the samples
tot_seqs_rrna_filt <- as.numeric(gen_stats_rrna_filt[,grep("total_sequences", colnames(gen_stats_rrna_filt))])
names(tot_seqs_rrna_filt) <- gsub("_fwd", "", gen_stats_rrna_filt$Sample)

# parse names
names(tot_seqs_rrna_filt)  <- paste("P",unlist(lapply(names(tot_seqs_rrna_filt), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## kegg read-based
setwd(paste(project_root,"/metatranscriptomics/kegg_diamond/combined_hits", sep = ""))

all_sample_info_files <- list.files()
total_kegg_counts <- numeric(length(all_sample_info_files))
total_ko_counts <- numeric(length(all_sample_info_files))

for(f in 1:length(all_sample_info_files)){
  # read in the mate-pair combined kegg hits 
  kegg_combined_hits <- fread(file   = all_sample_info_files[f], sep    = "\t", header = TRUE, quote  = "", fill   = FALSE, data.table = FALSE)
  
  # the number of rows is the total number of read-pairs (either R1, R2 or both) with a valid kegg gene hit
  total_kegg_counts[f] <- nrow(kegg_combined_hits)
  
  # number of read-pairs with a KO annotation, used for the KO analysis
  total_ko_counts[f] <- length(which(kegg_combined_hits$ko!=""))
}
names(total_kegg_counts) <- names(total_ko_counts) <- gsub(".txt", "", all_sample_info_files)

# parse names
names(total_kegg_counts)  <- paste("P",unlist(lapply(names(total_kegg_counts), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
names(total_ko_counts)  <- paste("P",unlist(lapply(names(total_ko_counts), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## metmarkdb read-based
setwd(paste(project_root,"/metatranscriptomics/metmarkdb_diamond/combined_hits", sep = ""))

all_sample_info_files <- list.files()
total_metmarkdb_counts <- numeric(length(all_sample_info_files))

for(f in 1:length(all_sample_info_files)){
  # read in the mate-pair combined metmarkdb hits 
  metmarkdb_combined_hits <- fread(file   = all_sample_info_files[f], sep    = "\t", header = FALSE, quote  = "", fill   = FALSE, data.table = FALSE)
  
  # the number of rows is the total number of read-pairs (either R1, R2 or both) with a valid kegg gene hit
  total_metmarkdb_counts[f] <- nrow(metmarkdb_combined_hits)
  
  # all of these hits have a marker gene annotation
}
names(total_metmarkdb_counts) <- gsub(".txt", "", all_sample_info_files)

# parse names
names(total_metmarkdb_counts)  <- paste("P",unlist(lapply(names(total_metmarkdb_counts), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

## phyloflash, read-based valid 16s counts in phyloflash output
# use the aggregated phyloflash output counts for this
setwd(paste(project_root,"/metatranscriptomics/phyloflash/downstream", sep = ""))
load("Otu_Tax_Tables_Parsed.RData")

phyloflash_total_counts <- colSums(otu_table_all, na.rm = TRUE)

# organize all stats similarly
tot_seqs <- tot_seqs[names(phyloflash_total_counts)]
tot_seqs_trimmed <- tot_seqs_trimmed[names(phyloflash_total_counts)]
tot_seqs_rrna_filt <- tot_seqs_rrna_filt[names(phyloflash_total_counts)]
total_kegg_counts <- total_kegg_counts[names(phyloflash_total_counts)]
total_ko_counts <- total_ko_counts[names(phyloflash_total_counts)]
total_metmarkdb_counts <- total_metmarkdb_counts[names(phyloflash_total_counts)]

# put into a matrix
mt_read_use_stats <- rbind(tot_seqs, tot_seqs_trimmed, tot_seqs_rrna_filt, total_kegg_counts, total_ko_counts, total_metmarkdb_counts, phyloflash_total_counts)
colnames(mt_read_use_stats) <- names(phyloflash_total_counts)
rownames(mt_read_use_stats) <- c("Total_read_pairs", "Total_read_pairs_trimmed", "Total_read_pairs_rrna_filtered", "Total_read_kegg_hits", "Total_read_kegg_hits_with_KO", "Total_read_metmarkdb_hits", "Total_read_phyloflash_counts")
mt_read_use_stats <- data.frame(mt_read_use_stats, stringsAsFactors = F, check.names = F)

# make additional tables as proportions
# metagenomics
mg_read_use_stats_prop <- mg_read_use_stats

# proportions to original number of read pairs
for(i in seq(from=2, to=nrow(mg_read_use_stats_prop))){
  mg_read_use_stats_prop[i,] <- as.numeric(mg_read_use_stats_prop[i,]) / as.numeric(mg_read_use_stats_prop[1,])
}

# proportions to trimmed reads
mg_read_use_stats_prop_trimmed <- mg_read_use_stats

# proportions to trimmed number of read pairs
for(i in seq(from=3, to=nrow(mg_read_use_stats_prop_trimmed))){
  mg_read_use_stats_prop_trimmed[i,] <- as.numeric(mg_read_use_stats_prop_trimmed[i,]) / as.numeric(mg_read_use_stats_prop_trimmed[2,])
}

# metatranscriptomics
mt_read_use_stats_prop <- mt_read_use_stats

# proportions to original number of read pairs
for(i in seq(from=2, to=nrow(mt_read_use_stats_prop))){
  mt_read_use_stats_prop[i,] <- as.numeric(mt_read_use_stats_prop[i,]) / as.numeric(mt_read_use_stats_prop[1,])
}

# proportions to trimmed reads
mt_read_use_stats_prop_trimmed <- mt_read_use_stats

# proportions to original number of read pairs
for(i in seq(from=3, to=nrow(mt_read_use_stats_prop_trimmed))){
  mt_read_use_stats_prop_trimmed[i,] <- as.numeric(mt_read_use_stats_prop_trimmed[i,]) / as.numeric(mt_read_use_stats_prop_trimmed[2,])
}

# proportions to rrna filtered reads reads
mt_read_use_stats_prop_rrna_filt <- mt_read_use_stats

# proportions to original number of read pairs
for(i in seq(from=4, to=nrow(mt_read_use_stats_prop_rrna_filt))){
  mt_read_use_stats_prop_rrna_filt[i,] <- as.numeric(mt_read_use_stats_prop_rrna_filt[i,]) / as.numeric(mt_read_use_stats_prop_rrna_filt[3,])
}

# change directory into the manuscript supplementary directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

# save the results tables into excel
# create a workbook
wb <- createWorkbook()

# metagenomics

# add worksheet
addWorksheet(wb, "MG_read_mapping")

# change the orientation
mg_table <- t(mg_read_use_stats)

# write the excel
writeData(wb, sheet = "MG_read_mapping", x = mg_table, startRow = 1, startCol = 1, rowNames = TRUE, colNames = TRUE)

# proportional table
mg_table_prop <- mg_read_use_stats_prop

# for all the functional and taxonomic annotations, use the proportions to the trimmed reads
mg_table_prop[3:nrow(mg_table_prop),] <- mg_read_use_stats_prop_trimmed[3:nrow(mg_read_use_stats_prop_trimmed),]

# change the orientation
mg_table_prop <- t(mg_table_prop)

# add worksheet
addWorksheet(wb, "MG_read_mapping_proportional")
# write the excel
writeData(wb, sheet = "MG_read_mapping_proportional", x = mg_table_prop, startRow = 1, startCol = 1, rowNames = TRUE, colNames = TRUE)

# metatranscriptomics
# change the orientation
mt_table <- t(mt_read_use_stats)

# add worksheet
addWorksheet(wb, "MT_read_mapping")
# write the excel
writeData(wb, sheet = "MT_read_mapping", x = mt_table, startRow = 1, startCol = 1, rowNames = TRUE, colNames = TRUE)

# proportional table
mt_table_prop <- mt_read_use_stats_prop

# for rrna filtered reads and taxonomy, use the proportions to trimmed reads
mt_table_prop[3,] <- mt_read_use_stats_prop_trimmed[3,]
mt_table_prop[7,] <- mt_read_use_stats_prop_trimmed[7,]

# for functional annotations, use the proportion to the rrna filtered reads
mt_table_prop[4:6,] <- mt_read_use_stats_prop_rrna_filt[4:6,]

# change the orientation
mt_table_prop <- t(mt_table_prop)

# add worksheet
addWorksheet(wb, "MT_read_mapping_proportional")
# write the excel
writeData(wb, sheet = "MT_read_mapping_proportional", x = mt_table_prop, startRow = 1, startCol = 1, rowNames = TRUE, colNames = TRUE)

# save the compiled excel workbook
saveWorkbook(wb, "MG_MT_read_mapping_stats.xlsx", overwrite = TRUE)

ses_info <- sessionInfo()

# save the read use stats also as an RData file for possible later use
save(mg_read_use_stats, mg_read_use_stats_prop, mg_read_use_stats_prop_trimmed, mt_read_use_stats, mt_read_use_stats_prop,
     mt_read_use_stats_prop_trimmed, mt_read_use_stats_prop_rrna_filt, mg_table, mg_table_prop, mt_table, mt_table_prop,ses_info, file = "mg_mt_read_use_stats.RData")
