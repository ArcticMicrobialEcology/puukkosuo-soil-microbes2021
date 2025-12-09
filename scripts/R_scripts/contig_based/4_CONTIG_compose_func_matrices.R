args <- commandArgs(trailingOnly = TRUE)

# input 1 needs to be the RData file for the kegg read based metagenomics normalized datas 
# input 2 needs to be the RData file for the kegg read based metranscriptomics normalized datas
# input 3 needs to be the downstream directory with the contig gene and transcript counts & annotations processed with a previous script

# e.g. 
# args <- character(3)
# args[1] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/metagenomics/kegg_diamond/Matrices_For_Downstream.RData"
# args[2] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/metatranscriptomics/kegg_diamond/Matrices_For_Downstream.RData"
# args[3] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/manu_rerun/downstream/contig_based"

# a small script to aggregate gene&transcript counts to the metabolic marker genes 
keep_items <- ls()
keep_items <- c(keep_items, "keep_items")

# load read-based filtered datas to get the sample metadata and all sample info
load(as.character(args[1]))
all_sample_info_mg <- all_sample_info

load(as.character(args[2]))
all_sample_info_mt <- all_sample_info

keep_items <- c(keep_items, "metadata", "all_sample_info_mg", "all_sample_info_mt")

# clean workspace, keep only the relevant objects
del_items <- ls()
del_items <- del_items[-which(del_items%in%keep_items)]
rm(list = del_items)

# get working directory
wd <- as.character(args[3])

# process metagenomics and metatranscriptomcis similarly
keep_items1 <- ls()
keep_items1 <- c(keep_items1, "keep_items1")
for(om in 1:2){
  setwd(wd)
  
  if(om==1){
    # load the compiled anvio mg results
    load("gene_counts_annotations.RData")
    
    # use the correct all sample info - metagenomics is processed first
    all_sample_info <- all_sample_info_mg
    
  }else{
    load("transcript_counts_annotations.RData")
    gene_annotations <- transcript_annotations
    gene_sample_counts <- transcript_sample_counts
    gene_info <- transcript_info
    
    # use the correct all sample info - metagenomics is processed first
    all_sample_info <- all_sample_info_mt
  }
  
  # multiple different types of annotations can be processed in the same way
  # currently processing only MetmarkDB - others may need to be parsed differently
  
  # parse the MetMarkDB genes
  metmark_db_genes <- strsplit(x = gene_annotations$metmarkdb_function, split = ";")
  metmark_db_genes <- unlist(lapply(metmark_db_genes, function(x) x[1]))
  metmark_db_genes[is.na(metmark_db_genes)]=""
  gene_annotations$metmarkdb_genes <- metmark_db_genes 
  
  # what sort of annotations we want to process - pick the columns
  process_cols <- which(colnames(gene_annotations)%in%c("metmarkdb_genes"))
  
  # store the results in lists
  count_datas <- list()
  rpk_datas <- list()
  rpm_datas <- list()
  tpm_datas <- list()
  
  count_sample_scale_factors <- list()
  rpk_sample_scale_factors <- list()
  
  keep_items2 <- ls()
  keep_items2 <- c(keep_items2, "keep_items2")
  for(j in 1:length(process_cols)){
    
    print(paste("Processing", colnames(gene_annotations)[process_cols[j]], "annotations"))
    
    # get all unique features in the data
    all_feat <- unique(unlist(strsplit(x = gene_annotations[,process_cols[j]], split = "; ")))
    if(any(all_feat%in%"")){all_feat <- all_feat[-which(all_feat%in%"")]}
    
    # summarize into feature level
    # RPK
    samples_rpk <- data.frame(matrix(nrow = length(all_feat), ncol = ncol(gene_sample_counts)))
    rownames(samples_rpk) <- all_feat
    colnames(samples_rpk) <- colnames(gene_sample_counts)
    
    # count
    samples_count <- samples_rpk
    
    # get gene lengths for RPK
    # anvio follows the convention of string indexing and splicing that is identical to the way one does it in Python or C. 
    # This means that the index of the first nucleotide in any contig should be 0. In other words, for a gene call that starts at the position xth position and ends at position yth position,
    # we start counting from x-1, and not from x (but we still end at y). The start and stop positions in the input file should comply with this criterion. 
    # but gene length in nucleotides should be stop-start now.
    gene_lengths <- gene_annotations$stop-gene_annotations$start
    
    # in kilobases
    gene_lengths_kb <- gene_lengths / 1000
    
    for(i in 1:nrow(samples_rpk)){
      locs_int <- grep(paste("\\b",rownames(samples_rpk)[i],"\\b", sep = ""), gene_annotations[,process_cols[j]])
      
      if(length(locs_int)>0){
        samples_count[i,] <- colSums(gene_sample_counts[locs_int,])
        samples_rpk[i,] <- colSums(gene_sample_counts[locs_int,] / gene_lengths_kb[locs_int])
      }else{
        samples_count[i,] <- rep(0, ncol(gene_sample_counts))
        samples_rpk[i,] <- rep(0, ncol(gene_sample_counts))
      }
    }
    
    
    # get scaling factors for rpm and tpm
    all_sample_info <- all_sample_info[colnames(samples_rpk)]
    rpk_scaling_factors <- numeric(ncol(samples_rpk))
    count_scaling_factors <- numeric(ncol(samples_count))
    
    for(i in 1:length(all_sample_info)){
      temp <- all_sample_info[[i]]
      rownames(temp) <- as.character(temp[,1])
      
      # get the total counts for each sample
      tot_counts <- as.numeric(temp["Total_Counts",2])
      # calculate the counts per million scaling factor
      count_scaling_factors[i] <- tot_counts / 10^6
      
      # get the total RPK for each sample
      tot_rpk <- as.numeric(temp["Total_RPK",2])
      # calculate the RPK per million scaling factor
      rpk_scaling_factors[i] <- tot_rpk / 10^6
    }
    
    names(count_scaling_factors) <- names(rpk_scaling_factors) <- names(all_sample_info)
    
    # to get rpm, divide the count for each gene with the count per million scaling factor for each sample
    samples_rpm <- samples_count
    for(i in 1:ncol(samples_rpm)){samples_rpm[,i] <- samples_rpm[,i]/count_scaling_factors[i]}
    
    # to get tpm, divide the rpk for each gene with the rpk per million scaling factor for each sample
    samples_tpm <- samples_rpk
    for(i in 1:ncol(samples_tpm)){samples_tpm[,i] <- samples_tpm[,i]/rpk_scaling_factors[i]}
    
    count_datas[[j]] <- samples_count
    rpk_datas[[j]] <- samples_rpk
    rpm_datas[[j]] <- samples_rpm
    tpm_datas[[j]] <- samples_tpm
    
    count_sample_scale_factors[[j]] <- count_scaling_factors
    rpk_sample_scale_factors[[j]] <- rpk_scaling_factors
    
    del_items <- ls()
    del_items <- del_items[-which(del_items%in%keep_items2)]
    rm(list = del_items)
  }
  
  # add suitable names
  namn <- "metmarkdb_diamond"
  names(tpm_datas) <- names(rpk_datas) <- names(rpm_datas) <-  names(count_datas) <- names(count_sample_scale_factors) <- names(rpk_sample_scale_factors) <- namn
  
  setwd(wd)
  if(om==1){ # metagenomics
    dir.create("metagenomics")
    setwd("metagenomics")
  } else { # metatranscriptomics
    dir.create("metatranscriptomics")
    setwd("metatranscriptomics")
  }
  
  # name the matrices
  mat_names <- character(length(count_datas))
  mat_names[grep("metmarkdb", namn)] <- "hits_samples"
  
  # save everything in suitable format for the downstream scripts
  for(i in 1:length(count_datas)){
    
    dir.create(names(count_datas)[i])
    setwd(names(count_datas)[i])
    assign(x = paste(mat_names[i], "_rpk", sep = ""), value = rpk_datas[[i]])
    assign(x = paste(mat_names[i], "_rpm", sep = ""), value = rpm_datas[[i]])
    assign(x = paste(mat_names[i], "_tpm", sep = ""), value = tpm_datas[[i]])
    assign(x = paste(mat_names[i], "_count", sep = ""), value = count_datas[[i]])
    
    assign(x = "count_scaling_factors", value = count_sample_scale_factors[[i]])
    assign(x = "rpk_scaling_factors", value = rpk_sample_scale_factors[[i]])
    
    write.table(x = get(paste(mat_names[i], "_rpk", sep = "")), file = "RPK_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    write.table(x = get(paste(mat_names[i], "_rpm", sep = "")), file = "RPM_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    write.table(x = get(paste(mat_names[i], "_tpm", sep = "")), file = "TPM_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    write.table(x = get(paste(mat_names[i], "_count", sep = "")), file = "Counts_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    
    save(list=c(paste(mat_names[i], "_rpk", sep = ""), paste(mat_names[i], "_rpm", sep = ""), paste(mat_names[i], "_tpm", sep = ""), 
                paste(mat_names[i], "_count", sep = ""), "all_sample_info", "count_scaling_factors", "rpk_scaling_factors"), file = "All_The_Output_Matrices_For_Downstream.RData")
    setwd("../")
  }
  
  del_items <- ls()
  del_items <- del_items[-which(del_items%in%keep_items1)]
  rm(list = del_items)
}

# print out session info
print("SessionInfo:")
sessionInfo()
