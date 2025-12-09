args <- commandArgs(trailingOnly = TRUE)

# the first input parameter needs to be the prodigal gene directory for the contigs
# the second input needs to be the the contig analysis main directory
# the third input parameter needs to be the anvio_project_name for the analysis

# e.g. 
# args <- character(3)
# args[1] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/manu_rerun/metagenomics/contig_based/prodigal_genes"
# args[2] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/manu_rerun/downstream/contig_based"
# args[3] <- "Oulanka_contig_based"

# a small script to parse all the functional and taxonomic annotations for the contig genes, get gene and transcript counts and link the annotations to the counts

# load libraries
library(grr)

# go to the prodigal gene directory 
setwd(as.character(args[1]))

# read-in all the identified prodigal genes in the contigs
all_genes <- read.csv("prodigal_genes_no_aaseq.txt", sep = "\t", header = T, quote = "", fill = F, stringsAsFactors = F)

# remove unnecessary columns here (source and version)
all_genes <- all_genes[, -c(which(colnames(all_genes)%in%c("source", "version")))]

# start adding more information after the last column
start_col <- ncol(all_genes)+1

# go through all the desired functional annotations - here could be several different
func_annot <- c("metmarkdb")
func_names <- character(0)

keep_items <- ls()
keep_items <- c(keep_items, "keep_items")

# read-in, match and save the functional annotations in the same data frame
for(i in 1:length(func_annot)){
  
  # read in the annotations
  filen <- paste("exported_functions_", func_annot[i], ".txt", sep = "")
  annotations_of_interest <- read.csv(file = filen, sep = "\t", header = T, quote = "", fill = F, stringsAsFactors = F)
  
  # get all matches based on gene callers id (prodigal)
  locs_matches <- matches(x = all_genes$gene_callers_id, y = annotations_of_interest$gene_callers_id, list = TRUE)
  names(locs_matches) <- all_genes$gene_callers_id
  
  # which genes don't have annotations
  no_hits <- which(lengths(locs_matches)==0)
  
  # which genes exactly one annotations
  unambiguous_hits <- which(lengths(locs_matches)==1)
  
  # which genes have multiple annotations
  multi_hits <- which(lengths(locs_matches)>1)
  
  # create a temporary variable for accession
  # parse in the annotations for the genes
  var <- rep("", nrow(all_genes))
  names(var) <- all_genes$gene_callers_id
  var[names(no_hits)] <- ""  
  var[names(unambiguous_hits)] <- annotations_of_interest$accession[unlist(locs_matches[names(unambiguous_hits)])]
  
  # process genes with multiple annotations separately
  if(length(multi_hits)>0){
    multi_list <- locs_matches[names(multi_hits)]
    multi_accession <- unlist(lapply(multi_list, function(x) paste(annotations_of_interest$accession[x], collapse = "; ")))
    var[names(multi_hits)] <- multi_accession
  }
  
  # add to the original gene data frame
  all_genes <-  cbind(all_genes, var)
  func_names <- c(func_names, paste(func_annot[i], "_accession", sep = ""))
  
  # create a temporary variable for additional functional annotations
  var <- rep("", nrow(all_genes))
  names(var) <- all_genes$gene_callers_id
  var[names(no_hits)] <- ""  
  var[names(unambiguous_hits)] <- annotations_of_interest$function.[unlist(locs_matches[names(unambiguous_hits)])]
  
  # process genes with multiple annotations separately
  if(length(multi_hits)>0){
    multi_accession <- unlist(lapply(multi_list, function(x) paste(annotations_of_interest$function.[x], collapse = "; ")))
    var[names(multi_hits)] <- multi_accession
  }
  
  # add to the original gene data frame
  all_genes <-  cbind(all_genes, var)
  func_names <- c(func_names, paste(func_annot[i], "_function", sep = ""))
  
  # clean workspace
  del_items <- ls()
  del_items <- del_items[-which(del_items%in%keep_items)]
  rm(list = del_items)
  
  # gc()
}

# parse column names
colnames(all_genes)[c(start_col:ncol(all_genes))] <- func_names

# save
setwd(as.character(args[2]))
write.table(x = all_genes, file = "master_annotation_functional_taxonomy.txt", quote = F, sep = "\t", col.names = T)

# read-in and parse the gene counts
setwd(as.character(args[1]))
gene_counts <- read.csv("gene_counts_all_samples.txt", sep = "\t", header = T, fill = F, stringsAsFactors = F, check.names = F, skip = 1)

# parse gene identifiers
anvio_project_name <- as.character(args[3])

# remove the project name from the gene identifier for matching
gene_counts$Geneid <- gsub(paste(anvio_project_name, "___", sep = ""), "", gene_counts$Geneid)

# separate the sample counts from the rest of the data and parse row and column names
gene_sample_counts <- gene_counts
gene_sample_counts <- gene_sample_counts[,-which(colnames(gene_sample_counts)%in%c("Geneid","Chr","Start","End","Strand","Length")),]
rownames(gene_sample_counts) <- gene_counts$Geneid

# parse sample names
colnames(gene_sample_counts) <- gsub("bam_files/", "", colnames(gene_sample_counts))
colnames(gene_sample_counts) <- gsub(".bam", "", colnames(gene_sample_counts))

# save also the additional data related to the mapping
gene_info <- gene_counts[,which(colnames(gene_counts)%in%c("Geneid","Chr","Start","End","Strand","Length")),]
rownames(gene_info) <- gene_info$Geneid
rm(gene_counts)# not needed anymore

# get annotations for the genes and match with the gene count data
rownames(all_genes) <- all_genes$gene_callers_id
gene_annotations <- all_genes[rownames(gene_sample_counts),]

# now here we have the gene counts for the samples and annotations for the same genes to be used 
# for downstream analysis

# save everything
setwd(as.character(args[2]))
write.table(x = gene_sample_counts, file = "gene_sample_counts.txt", quote = F, sep = "\t", col.names = T)
write.table(x = gene_info, file = "gene_info.txt", quote = F, sep = "\t", col.names = T)
write.table(x = gene_annotations, file = "gene_annotations.txt", quote = F, sep = "\t", col.names = T)
save(gene_sample_counts, gene_info, gene_annotations, file = "gene_counts_annotations.RData")

rm(gene_annotations)# not needed anymore
rm(gene_info)# not needed anymore
rm(gene_sample_counts)# not needed anymore

# process transcript counts similarly
# read-in the transcript counts
setwd(as.character(args[1]))
transcript_counts <- read.csv("transcript_counts_all_samples.txt", sep = "\t", header = T, fill = F, stringsAsFactors = F, check.names = F, skip = 1)

# remove the project name from the gene identifier for matching
transcript_counts$Geneid <- gsub(paste(anvio_project_name, "___", sep = ""), "", transcript_counts$Geneid)

# separate the samples from the rest of the data and parse row and column names
transcript_sample_counts <- transcript_counts
transcript_sample_counts <- transcript_sample_counts[,-which(colnames(transcript_sample_counts)%in%c("Geneid","Chr","Start","End","Strand","Length")),]
rownames(transcript_sample_counts) <- transcript_counts$Geneid

# parse sample names
colnames(transcript_sample_counts) <- gsub("bam_files_metat/", "", colnames(transcript_sample_counts))
colnames(transcript_sample_counts) <- gsub(".bam", "", colnames(transcript_sample_counts))

# save also the additional data related to the mapping
transcript_info <- transcript_counts[,which(colnames(transcript_counts)%in%c("Geneid","Chr","Start","End","Strand","Length")),]
rownames(transcript_info) <- transcript_info$Geneid
rm(transcript_counts) # not needed anymore

# get annotations
transcript_annotations <- all_genes[rownames(transcript_sample_counts),]

# now here we have the transcript counts for identified genes in the samples and annotations for the same genes to be used 
# for downstream analysis

# save everything
setwd(as.character(args[2]))
write.table(x = transcript_sample_counts, file = "transcript_sample_counts.txt", quote = F, sep = "\t", col.names = T)
write.table(x = transcript_info, file = "transcript_info.txt", quote = F, sep = "\t", col.names = T)
write.table(x = transcript_annotations, file = "transcript_annotations.txt", quote = F, sep = "\t", col.names = T)
save(transcript_sample_counts, transcript_info, transcript_annotations, file = "transcript_counts_annotations.RData")
# done

# print out session info
print("SessionInfo:")
sessionInfo()
