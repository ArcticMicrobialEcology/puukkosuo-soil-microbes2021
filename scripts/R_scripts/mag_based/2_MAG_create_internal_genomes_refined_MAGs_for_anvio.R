args <- commandArgs(trailingOnly = TRUE)

# the first input parameter needs to be the output directory where the internal_genomes.txt containing the MAGs 
# for the different coassemblies for anvio will be saved
# the rest of the input parameters / arguments need to be the bin_by_bin directories of the user refined Metabat2 MAGs for
# the different coassemblies produced by anvi-summarize

# a small script to produce an "internal genomes" file for the MAGs from the different coassemblies for anvi-dereplicate-genomes

# get all the refined MAGs
all_refined_bins <- list()
for(i in 2:length(args)){
  setwd(args[i])
  all_refined_bins[[i-1]] <- list.files()
}

# parse coassembly id, assumes the same directory structure always
# if there are changes to this practice, check this
# same with all directory parsing performed by this script
coassembly_names <- strsplit(x = args[2:length(args)], split = "/")
coassembly_names <- unlist(lapply(coassembly_names, function(x) tail(x,5)[1]))
names(all_refined_bins) <- coassembly_names

# required information for the internal genomes file for anvio:
# name bin_id	collection_id	profile_db_path	contigs_db_path
internal_genomes <- data.frame(matrix(nrow = sum(lengths(all_refined_bins)), ncol = 5))

# start processing
ind <- 1
for(i in 1:length(all_refined_bins)){
  
  # save the information for all the coassemblies
  internal_genomes[ind:(ind+length(all_refined_bins[[i]])-1),1] <- paste(names(all_refined_bins)[i], all_refined_bins[[i]], sep = "_")
  internal_genomes[ind:(ind+length(all_refined_bins[[i]])-1),2] <- all_refined_bins[[i]]
  
  # get the profile database paths for the MAGs
  # assumes a specific structure (see args)
  coassembly_dir <- unlist(strsplit(x = args[i+1], split = "/"))
  coassembly_dir <- rev(rev(coassembly_dir)[-c(1:4)])
  coassembly_dir <- c(coassembly_dir, "SAMPLES-MERGED", "PROFILE.db")
  coassembly_dir <- paste(coassembly_dir, collapse = "/")
  internal_genomes[ind:(ind+length(all_refined_bins[[i]])-1),4] <- rep(coassembly_dir, length(all_refined_bins[[i]]))
  
  # get the contig database paths for the MAGs
  # assumes a specific structure (see args)
  coassembly_dir <- unlist(strsplit(x = args[i+1], split = "/"))
  coassembly_dir <- rev(rev(coassembly_dir)[-c(1:4)])
  coassembly_dir <- paste(coassembly_dir, collapse = "/")
  coassembly_files <- list.files(path = coassembly_dir)
  contig_db_file <- coassembly_files[grep("_contigs.db",coassembly_files)]
  coassembly_dir <- paste(coassembly_dir, "/", contig_db_file, sep = "")
  internal_genomes[ind:(ind+length(all_refined_bins[[i]])-1),5] <- rep(coassembly_dir, length(all_refined_bins[[i]]))
  
  ind <- ind+length(all_refined_bins[[i]])
  #ind <- length(all_refined_bins[[i]])+1
}

# same collection_id for all mags
internal_genomes[,3] <- rep("metabat2_refined_mags", nrow(internal_genomes))
colnames(internal_genomes) <- c("name","bin_id","collection_id","profile_db_path","contigs_db_path")

# save the internal genomes list here
setwd(args[1])
write.table(x = internal_genomes, file = "internal_genomes_refined_MAGs.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
