args <- commandArgs(trailingOnly = TRUE)

# the first input parameter needs to be the directory for metagenomics with all the result subdirectories
# the second input parameter needs to be the directory for metatranscriptomics with all the result subdirectories
# the third input parameter needs to be the preprocessed metadata file

# e.g.
# args <- character(3)
# args[1] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/manu_rerun/downstream/contig_based/metagenomics"
# args[2] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/manu_rerun/downstream/contig_based/metatranscriptomics"
# args[3] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/manu_rerun/metadata/Study_Metadata.RData"

load(as.character(args[3]))
rownames(metadata) <- metadata$`Short code`

# a small script to filter the metabolic marker gene datas for lowly expressed features and to parse the 
# sample names and orders to match the metadata

for(om in 1:2){
  
  print(paste("Processing data from directory:", as.character(args[om])))
  
  # change working directory
  setwd(as.character(args[om]))
  
  # load and process the metmarkdb data
  setwd("metmarkdb_diamond")
  load("All_The_Output_Matrices_For_Downstream.RData")
  
  # parse sample names and order to match metadata
  colnames(hits_samples_tpm) <- paste("P",unlist(lapply(colnames(hits_samples_tpm), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
  hits_samples_tpm <- hits_samples_tpm[,rownames(metadata)]
  
  colnames(hits_samples_count) <- paste("P",unlist(lapply(colnames(hits_samples_count), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
  hits_samples_count <- hits_samples_count[,rownames(metadata)]
  
  colnames(hits_samples_rpm) <- paste("P",unlist(lapply(colnames(hits_samples_rpm), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
  hits_samples_rpm <- hits_samples_rpm[,rownames(metadata)]
  
  names(rpk_scaling_factors) <- paste("P",unlist(lapply(names(rpk_scaling_factors), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
  rpk_scaling_factors <- rpk_scaling_factors[rownames(metadata)]
  
  names(count_scaling_factors) <- paste("P",unlist(lapply(names(count_scaling_factors), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
  count_scaling_factors <- count_scaling_factors[rownames(metadata)]
  
  # filter similarly to the KEGG data
  
  # some expression is defined as count 5
  min_count <- 5
  
  # minimum amount of sample is defined as 6, corresponding to the amount of samples in the vegetation
  # cluster or the smallest subgroup of grazing*snowtreatment
  min_samples <- 6
  
  # filter tpm data and rpm data based on rpm thres
  loc_min_sample <- which.min(count_scaling_factors)
  locs_limit_vals <- which(hits_samples_count[,loc_min_sample]==min_count)
  if(length(locs_limit_vals)<1){
    cand_vals <-  hits_samples_count[,loc_min_sample]
    
    # don't accept lower counts
    cand_vals[which(cand_vals<min_count)] <- NA
    
    # pick a new threshold
    locs_limit_vals <- which.min(abs(cand_vals - min_count))
  }
  
  # use the rpm data for filtering, where the counts for features have been adjuted to the library size
  thres_expr <- unique(hits_samples_rpm[locs_limit_vals, loc_min_sample])
  keep <- rowSums(hits_samples_rpm > thres_expr ) >= min_samples
  
  metmark_tpm_data <- hits_samples_tpm[keep,]
  metmark_rpm_data <- hits_samples_rpm[keep,]
  metmark_count_data <- hits_samples_count[keep,]
  
  # save everything
  save(metmark_tpm_data, metmark_rpm_data, metmark_count_data, hits_samples_tpm, hits_samples_rpm, hits_samples_rpk, hits_samples_count, all_sample_info, count_scaling_factors, rpk_scaling_factors,
       metadata, min_count, min_samples, thres_expr, file = "Matrices_For_Downstream.RData")
  
}

# print out session info
print("SessionInfo:")
sessionInfo()
