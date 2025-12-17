args <- commandArgs(trailingOnly = TRUE)

# the first parameter needs to be the prepared (by the anvio contig workflow) reformatting report from the concatenated contig assembly directory, with separate conditions named
# the second parameter needs to be directory for the prodigal gene counts for the concatenated contigs processing
# the third parameter needs to be the output directory for the mapping files
# the rest of the parameters starting from four need to be the prodigal gene directories for the separate coassemblies (one for each coassembly) for mag-based anvio workflow

# a small script to map the gene and contig annotations from the processing of separate coassemblies to the 
# the processing of the concatenated coassembly (where all contigs have been concatenated into one database)
# to allow for the use of gene and transcript counts from the concatenated database for the MAGs

# load libraries
library(grr)

# load the reformattings report
reformat_report <- read.csv(args[1], sep = "\t", header = F)
colnames(reformat_report) <- c("new","old") # should be the other way around according to anvio documentation, but is not. A bug?

# parse condition in report
reformat_report_condition <- unlist(lapply(strsplit(x = as.character(reformat_report$old), split = "_"), function(x) x[3]))
reformat_report$condition <- reformat_report_condition
reformat_report$old <- unlist(lapply(strsplit(x = as.character(reformat_report$old), split = "_"), function(x) paste(x[1:2], collapse = "_")))

# load the concantenated results
setwd(args[2])
load("gene_counts_annotations.RData")

# genes
gene_sample_counts_concatenated <- gene_sample_counts
gene_info_concatenated <- gene_info
gene_annotations_concatenated <- gene_annotations

# load the results from the different coassemblies and map them to the coassembly
# put the mapping files in this list
mapping_files <- list()

# get coassembly names - assumes a specific directory structure - needs to be corrected if modified
coassembly_names <- strsplit(x = args[4:length(args)], split = "/")
coassembly_names <- unlist(lapply(coassembly_names, function(x) tail(x,2)[1]))

for(i in 4:length(args)){
  
  print(paste("Processing condition / coassembly:", coassembly_names[i-3]))
  setwd(args[i])
  
  # read-in the gene counts
  gene_counts <- read.csv("gene_counts_all_samples.txt", sep = "\t", header = T, fill = F, stringsAsFactors = F, check.names = F, skip = 1)
  
  # parse gene identifiers
  anvio_project_name <- paste("Oulanka_", coassembly_names[i-3], sep = "")
  
  # remove the project name from the gene identifier for matching
  gene_counts$Geneid <- gsub(paste(anvio_project_name, "___", sep = ""), "", gene_counts$Geneid)
  
  # separate the sample counts from the rest of the data and parse row and column names
  gene_sample_counts <- gene_counts
  gene_sample_counts <- gene_sample_counts[,-which(colnames(gene_sample_counts)%in%c("Geneid","Chr","Start","End","Strand","Length")),]
  rownames(gene_sample_counts) <- gene_counts$Geneid
  colnames(gene_sample_counts) <- gsub("bam_files/", "", colnames(gene_sample_counts))
  colnames(gene_sample_counts) <- gsub(".bam", "", colnames(gene_sample_counts))
  
  gene_info <- gene_counts[,which(colnames(gene_counts)%in%c("Geneid","Chr","Start","End","Strand","Length")),]
  rownames(gene_info) <- gene_info$Geneid
  gene_info$Chr <- gsub(paste("_",coassembly_names[i-3], sep = ""),"", gene_info$Chr)
  rm(gene_counts)# not needed anymore
  
  # get the relevant contigs to the coassembly from the reformatting report
  reformat_report_coassembly <- reformat_report[which(reformat_report$condition%in%coassembly_names[i-3]),]
  
  # get the mapping information
  locs_contig_matches <- matches(x = gene_info$Chr, y = reformat_report_coassembly$old , list = TRUE)
  
  # put-in the new names
  names(locs_contig_matches) <- gene_info$Geneid
  gene_info$New_Chr <- character(nrow(gene_info))
  gene_info[names(locs_contig_matches),which(colnames(gene_info)%in%"New_Chr")] <- reformat_report_coassembly$new[unlist(locs_contig_matches)]
  
  # new gene names
  conc_contig_coassembly <- gene_info_concatenated[which(gene_info_concatenated$Chr%in%gene_info$New_Chr),]
  conc_contig_coassembly <- conc_contig_coassembly[order(conc_contig_coassembly$Chr, conc_contig_coassembly$Length),]
  gene_info <- gene_info[order(gene_info$New_Chr, gene_info$Length),]
  
  # do some checks
  check1 <- all(conc_contig_coassembly$Chr==gene_info$New_Chr)
  check2 <- all(conc_contig_coassembly$Length==gene_info$Length)
  print(paste("All contig IDs match:", check1))
  print(paste("All gene lengths match:", check2))
  
  # put the new gene names
  gene_info$New_Geneid <- conc_contig_coassembly$Geneid
  
  # save the mapping file
  save(gene_info, file = "Gene_Contig_Mapping_To_Conc_Coassem.RData")
  mapping_files[[i-3]] <- gene_info
}
names(mapping_files) <- coassembly_names

# save the mapping files
setwd(as.character(args[3]))
save(mapping_files, coassembly_names, file = "Gene_Contig_Mapping_Sep_Coassem_Conc_Coassem.RData")

# print out session info
sessionInfo()
