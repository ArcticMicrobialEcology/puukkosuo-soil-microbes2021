args <- commandArgs(trailingOnly = TRUE)

# the first input parameter needs to be the anvio MAG dereplication result directory containing the dereplication report
# the second input parameter needs to be the output directory where the final dereplicated MAGs are copied to (needs to exist)
# the third input parameter needs to be the MAG base directory

# a script to gather all the dereplicated final mags and all the related information to them and copy them into a separate
# directory for easier downstream processing

# load libraries
library(dplyr)
library(tidyr)

# read-in the dereplication report
setwd(as.character(args[1]))
dereplication_report <- read.csv("CLUSTER_REPORT.txt", sep = "\t")

# store some extra information in these variables
bin_conditions <- character(nrow(dereplication_report))
bin_names <- character(nrow(dereplication_report))
bin_dirs <- character(nrow(dereplication_report))

# go throught each bin cluster in the report
for(bin in 1:nrow(dereplication_report)){
  
  # get the representative bin
  representative_bin <- dereplication_report$representative[bin]
  
  # get the condition for the representative bin and the name for the bin in the condition
  representative_bin_split <- strsplit(x = representative_bin, split = "_")
  cut_ind <- grep("metabat2", representative_bin_split[[1]])
  representative_bin_condition <- paste(representative_bin_split[[1]][1:(cut_ind-1)], collapse = "_")
  representative_bin_name <- paste(representative_bin_split[[1]][cut_ind:length(representative_bin_split[[1]])], collapse = "_")
  
  # the directory in which the representative bin should be located
  representative_bin_dir <- paste(as.character(args[3]), "/", representative_bin_condition, "/MAGs/metabat2/summary_refined_mags/bin_by_bin", sep = "")
  
  # check that such a bin exists
  setwd(representative_bin_dir)
  cond_files <- list.files()
  check1 <- representative_bin_name%in%cond_files
  
  print(paste("Processing MAG cluster:", dereplication_report$cluster[bin]))
  print(paste("Representative bin:", representative_bin))
  print(paste("From condition:", representative_bin_condition))
  print(paste("Bin name in condition:",representative_bin_name))
  print(paste("Found representative bin in condition directory:", check1))
  
  # create a directory for the bin cluster in the final MAG directory
  dir.create(paste(as.character(args[2]), "/", dereplication_report$cluster[bin], sep=""))
  
  # copy the representative bin files to the final MAG directory
  setwd(representative_bin_name)
  representative_bin_files <- list.files()
  file.copy(from = representative_bin_files, to = paste(as.character(args[2]), "/", dereplication_report$cluster[bin], sep=""))
  
  # save some extra information for the bin for easier downstream processing
  bin_conditions[bin] <- representative_bin_condition
  bin_names[bin] <- representative_bin_name
  bin_dirs[bin] <- gsub("/bin_by_bin", "", representative_bin_dir)
  
  print("*****************************************")
  print("*****************************************")
}

# store the extra variables in the dereplication report
dereplication_report$condition <- bin_conditions
dereplication_report$bin_name_in_condition <- bin_names
dereplication_report$base_dir <- bin_dirs

# make some annotations and quantative matrices for downstream analysis
# gtdbtk
unique_base_dirs <- unique(dereplication_report$base_dir)
gtdbtk_reports <- list()

# compile a single taxonomy matrix for all the clusters
gtdbtk_report <- data.frame(matrix(ncol = 7, nrow = nrow(dereplication_report)), stringsAsFactors = F)
rownames(gtdbtk_report) <- dereplication_report$cluster
colnames(gtdbtk_report) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for(dir in 1:length(unique_base_dirs)){
  
  print(paste("Processing GTDB-TK taxonomy from directory:", unique_base_dirs[dir]))
  
  setwd(unique_base_dirs[dir])
  setwd("gtdbtk")
  gtdbtk_temp <- read.csv("gtdbtk.bac120.summary.tsv", sep = "\t")
  gtdbtk_temp$bin_name <- gsub("-contigs", "", gtdbtk_temp$user_genome)
  
  # input taxonomy to the common taxonomy table for all final MAGs 
  locs_dir_included <- which(dereplication_report$base_dir%in%unique_base_dirs[dir])
  temp_bins <- dereplication_report$bin_name_in_condition[locs_dir_included]
  
  # get relevant gtbk-tk results
  rownames(gtdbtk_temp) <- gtdbtk_temp$bin_name
  gtdbtk_temp <- gtdbtk_temp[temp_bins,]
  
  # split taxonomy
  temp_taxonomy <- strsplit(x = gtdbtk_temp$classification, split = ";")
  
  # check lengths
  check1 <- all(lengths(temp_taxonomy)==7)
  print(paste("All MAGs have 7 levels of taxonomy:", check1))
  
  # convert into data frame
  temp_taxonomy <- do.call(rbind.data.frame, temp_taxonomy)
  print(paste("All Domain annotaions seem to be Domain annotations:", length(grep("d__", temp_taxonomy[,1]))==nrow(temp_taxonomy)))
  print(paste("All Phylum annotaions seem to be Phylum annotations:", length(grep("p__", temp_taxonomy[,2]))==nrow(temp_taxonomy)))
  print(paste("All Class annotaions seem to be Class annotations:", length(grep("c__", temp_taxonomy[,3]))==nrow(temp_taxonomy)))
  print(paste("All Order annotaions seem to be Order annotations:", length(grep("o__", temp_taxonomy[,4]))==nrow(temp_taxonomy)))
  print(paste("All Family annotaions seem to be Family annotations:", length(grep("f__", temp_taxonomy[,5]))==nrow(temp_taxonomy)))
  print(paste("All Genus annotaions seem to be Genus annotations:", length(grep("g__", temp_taxonomy[,6]))==nrow(temp_taxonomy)))
  print(paste("All Species annotaions seem to be Species annotations:", length(grep("s__", temp_taxonomy[,7]))==nrow(temp_taxonomy)))
  
  # parse annotations
  temp_taxonomy[,1] <- gsub("d__", "", temp_taxonomy[,1])
  temp_taxonomy[,2] <- gsub("p__", "", temp_taxonomy[,2])
  temp_taxonomy[,3] <- gsub("c__", "", temp_taxonomy[,3])
  temp_taxonomy[,4] <- gsub("o__", "", temp_taxonomy[,4])
  temp_taxonomy[,5] <- gsub("f__", "", temp_taxonomy[,5])
  temp_taxonomy[,6] <- gsub("g__", "", temp_taxonomy[,6])
  temp_taxonomy[,7] <- gsub("s__", "", temp_taxonomy[,7])
  colnames(temp_taxonomy) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(temp_taxonomy) <- gtdbtk_temp$bin_name
  
  gtdbtk_report[locs_dir_included,] <- temp_taxonomy
  gtdbtk_reports[[dir]] <- gtdbtk_temp
  
  print("*****************************************")
  print("*****************************************")
}
names(gtdbtk_reports) <- unique_base_dirs

# # get the number of copies for the marker genes for all MAG clusters
# load marker gene metadata
load("/scratch/project_2009164/2_OULANKA/Tommi/final/databases/metmarkdb/Gene_metadata.RData")

# put in a data frame
markergene_copies_clusters <- data.frame(matrix(nrow = nrow(gene_metadata), ncol = nrow(dereplication_report)))
rownames(markergene_copies_clusters) <- rownames(gene_metadata)
colnames(markergene_copies_clusters) <- dereplication_report$cluster

# get the copy numbers for each MAG cluster
all_bin_dirs <- paste(dereplication_report$base_dir, "/bin_by_bin/", dereplication_report$bin_name_in_condition, sep = "")
for(bin in 1:length(all_bin_dirs)){
  
  # change directory
  setwd(all_bin_dirs[bin])
  
  # read-in the gene file
  all_bin_files <- list.files()
  bin_gene_file <- read.csv(all_bin_files[grep("-gene_calls.txt", all_bin_files)], sep = "\t")
  
  # marker genes
  all_bin_marker_genes <- as.character(bin_gene_file$MetMarkDB)
  
  # get only the marker genes
  all_bin_marker_genes <- strsplit(x = all_bin_marker_genes, split = ";")
  all_bin_marker_genes <- unlist(lapply(all_bin_marker_genes, function(x) x[1]))
  all_bin_marker_genes <- as.character(na.omit(all_bin_marker_genes))
  if(any(all_bin_marker_genes == "")) {
    all_bin_marker_genes <- all_bin_marker_genes[-which(all_bin_marker_genes=="")]
  }
  markergene_copies_clusters[names(table(all_bin_marker_genes)), bin] <- table(all_bin_marker_genes)
}

# set NA as 0
markergene_copies_clusters[is.na(markergene_copies_clusters)] <- 0

# save everything
setwd(as.character(args[2]))

save(dereplication_report, gtdbtk_report, gtdbtk_reports, markergene_copies_clusters,
     file="Final_MAG_Taxonomy_Mark_Gene.RData")

# print - out session info
sessionInfo()
