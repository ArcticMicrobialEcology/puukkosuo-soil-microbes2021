args <- commandArgs(trailingOnly = TRUE)
# the first input parameter needs to be the "sample_details" directory for the samples resulting from the Rscript "2a_READ_combine_KEGG_alignment_paired_reads.R"
# the second parameter needs to be the directory for the KEGG alignment hits combined for the paired reads and aggregated to the gene level from the same script.
# the third parameter needs to be the result directory (needs to exist) where the KO summarized KEGG alignment result matrices are saved.

# e.g. 
# args <- character(3)
# args[1] <- "/scratch/project_2007998/AGROBIO_Tommi/metagenomics/kegg_diamond/sample_details"
# args[2] <- "/scratch/project_2007998/AGROBIO_Tommi/metagenomics/kegg_diamond/combined_hits_gene_aggregated"
# args[3] <- "/scratch/project_2007998/AGROBIO_Tommi/metagenomics/kegg_diamond"

# this is a small script to summarize the the KEGG alignment hits combined for the paired reads to KO level and transform into matrices for downstream processing.

# get sample specific total counts and total rpks, calculated based on KEGG alignments
setwd(as.character(args[1]))
all_sample_info_files <- list.files()
all_sample_info <- list()

for(f in 1:length(all_sample_info_files)){
  all_sample_info[[f]] <- read.csv(all_sample_info_files[f], sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F)
}
names(all_sample_info) <- gsub(".txt", "", all_sample_info_files)

# read in the the KEGG alignment hits combined for the paired reads and aggregated to gene level
setwd(args[2])
all_files <- list.files()

# read into a list
agg_hits_list <- list()
for(f in 1:length(all_files)){
  kegg_hits <- read.csv(all_files[f], sep = "\t", header = T, quote = "", fill = F, stringsAsFactors = F)
  
  # add sample information
  kegg_hits$sample <- rep(gsub(".txt", "", all_files[f]), nrow(kegg_hits))
  agg_hits_list[[f]] <- kegg_hits
  rm(kegg_hits)
}
names(agg_hits_list) <- gsub(".txt", "", all_files)

# combine hits from the different samples
print("Combining hits from the different samples")

# this will take a while
hits_all_samples <- do.call("rbind", agg_hits_list)

print("All samples combined")
print("Summarizing the detected hits to KO group level and converting into a matrix form for RPK and counts. Calculating TPM from the RPK.")

# get all unique KO groups in the data
all_pos_ko_data <- unique(unlist(strsplit(x = hits_all_samples$ko, split = ";")))
if(any(all_pos_ko_data=="")){all_pos_ko_data <- all_pos_ko_data[-which(all_pos_ko_data=="")]}

# summarize into KO group level
# RPK
ko_samples_rpk <- data.frame(matrix(nrow = length(all_pos_ko_data), ncol = length(agg_hits_list)))
rownames(ko_samples_rpk) <- all_pos_ko_data
colnames(ko_samples_rpk) <- names(agg_hits_list)

# count
ko_samples_count <- ko_samples_rpm <- ko_samples_rpk

# summarize hits for each sample and also record if a hit is calculated several times for the samples
rpk_cor_factors <- list()
count_cor_factors <- list()

for(i in 1:length(agg_hits_list)){
  temp <- agg_hits_list[[i]]
  
  for(j in 1:nrow(ko_samples_rpk)){
    locs_ko <- grep(paste("\\b",rownames(ko_samples_rpk)[j],"\\b", sep = ""), temp$ko)
    
    if(length(locs_ko)>0){
      temp2 <- temp[locs_ko,]
      
      ko_samples_rpk[j,i] <- sum(as.numeric(temp2$rpk))
      ko_samples_count[j,i] <- sum(as.numeric(temp2$count))
      ko_samples_rpm[j,i] <- sum(as.numeric(temp2$count))
    }
    rm(locs_ko)
    rm(temp2)
  }
  
  # as TPM = GENErpk (or FEATURErpk) / (SUM(SAMPLErpk) / 10^6)
  # Since (SUM(SAMPLErpk) / 10^6) has already been estimated from KEGG alignments 
  # (assuming all prokaryotic genes present in the sample are represented in the KEGG database),
  # protein hits assigned to multiple KO groups may have been counted multiple times.
  # To correct for this, count how many times each such protein was used per sample,
  # sum the RPK values of proteins counted more than once, 
  # and adjust the total SUM(SAMPLErpk) accordingly before TPM calculation.
  
  # calculate the additional RPK induced by multicounting to correct the RPK scaling factors for each sample
  # to retain comparability between sample TPM values even when multicounting the protein hits for several KO groups.
  
  # get the multimapped proteins
  locs_multi <- grep(";", temp$ko)
  
  # each of these gene/protein hits is counted more than once
  temp_multi <- temp[locs_multi,]
  
  # how many times - count the separators - this will directly result in the "extra counts"
  nr_multi <- unlist(lapply(temp_multi$ko, function(x) length(gregexpr(";", x)[[1]])))
  
  # additional RPK and counts from the multicounts for each sample
  rpk_sum_cor_factor <- sum(as.numeric(temp_multi$rpk) * nr_multi)
  count_sum_cor_factor <- sum(as.numeric(temp_multi$count) * nr_multi)
  
  rpk_cor_factors[[i]] <- rpk_sum_cor_factor
  count_cor_factors[[i]] <- count_sum_cor_factor
  
  rm(count_sum_cor_factor)
  rm(rpk_sum_cor_factor)
  rm(temp_multi)
  rm(nr_multi)
  rm(locs_multi)
  rm(temp)
}
names(rpk_cor_factors) <- names(agg_hits_list)
names(count_cor_factors) <- names(agg_hits_list)

# set NA to zeroes, change this if needed
ko_samples_rpk[is.na(ko_samples_rpk)] <- 0
ko_samples_rpm[is.na(ko_samples_rpm)] <- 0
ko_samples_count[is.na(ko_samples_count)] <- 0

# now, transform the RPK to TPM by using the RPK scaling factors for the samples from the KEGG alignment. 
# the KEGG alignment is a fairly good estimate of all the prokaryotic genes in the samples and thus the total gene counts and total RPK in the samples.

# also correct scaling factors with the sum RPK of multicounted gene/protein hits for each sample
all_sample_info <- all_sample_info[colnames(ko_samples_rpk)]
rpk_scaling_factors <- numeric(ncol(ko_samples_rpk))
count_scaling_factors <- numeric(ncol(ko_samples_rpm))

# produce also RPM data similarly to RPK data
for(i in 1:length(all_sample_info)){
  temp <- all_sample_info[[i]]
  rownames(temp) <- as.character(temp[,1])
  
  # get the total RPK for each sample
  tot_rpk <- as.numeric(temp["Total_RPK",2])
  # get the total counts for each sample
  tot_counts <- as.numeric(temp["Total_Counts",2])
  
  # correct with the additional RPK from the multicounted hits
  tot_rpk <- tot_rpk + as.numeric(rpk_cor_factors[[i]])
  # correct with the additional counts from the multicounted hits
  tot_counts <- tot_counts + as.numeric(count_cor_factors[[i]])
  
  # calculate the RPK per million scaling factor
  rpk_scaling_factors[i] <- tot_rpk / 10^6
  # calculate the counts per million scaling factor
  count_scaling_factors[i] <- tot_counts / 10^6
    
}
names(rpk_scaling_factors) <- colnames(ko_samples_rpk)
names(count_scaling_factors) <- colnames(ko_samples_rpm)

# to get TPM, divide the RPK for each gene with the RPK per million scaling factor for each sample
ko_samples_tpm <- ko_samples_rpk
for(i in 1:ncol(ko_samples_tpm)){ko_samples_tpm[,i] <- ko_samples_tpm[,i]/rpk_scaling_factors[i]}

# to get RPM, divide the counts for each gene with the per million scaling factor for each sample
for(i in 1:ncol(ko_samples_rpm)){ko_samples_rpm[,i] <- ko_samples_rpm[,i]/count_scaling_factors[i]}

# these matrices can be taken forward for analysis.
print("Saving final matrices")

setwd(args[3])
write.table(x = hits_all_samples, file = "Combined_Hits_All_Samples.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = ko_samples_rpk, file = "RPK_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(x = ko_samples_rpm, file = "RPM_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(x = ko_samples_tpm, file = "TPM_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(x = ko_samples_count, file = "Counts_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

save(ko_samples_rpk, ko_samples_rpm, ko_samples_tpm, ko_samples_count, all_sample_info, rpk_cor_factors, count_cor_factors, rpk_scaling_factors, count_scaling_factors, file = "All_The_Output_Matrices_For_Downstream.RData")

# print out session info
print("SessionInfo:")
sessionInfo()