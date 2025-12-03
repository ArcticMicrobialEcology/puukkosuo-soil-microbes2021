args <- commandArgs(trailingOnly = TRUE)
# the first input parameter needs to be the "sample_details" directory for the samples resulting from the Rscript "2a_READ_combine_KEGG_alignment_paired_reads.R"
# the second parameter needs to be the result directory for the DIAMOND MetMarkDB alignment results
# the third parameter needs to be the output directory (needs to exist) where the subdirectories for the results are created (if they don't exist) and saved
# the fourth parameter needs to be if the MetMarkDB hits should be filtered (1) or not (0) with the defined default limits.
# currently the default limits: e-value 10^-6 and bitscore 50, change if needed

# e.g.
# args <- character(4)
# args[1] <- "/scratch/project_2007998/AGROBIO_Tommi/metagenomics/kegg_diamond/sample_details"
# args[2] <- "/scratch/project_2007998/AGROBIO_Tommi/metagenomics/metmarkdb_diamond/raw"
# args[3] <- "/scratch/project_2007998/AGROBIO_Tommi/metagenomics/metmarkdb_diamond"
# args[4] <- "0"

# this is a small script to filter, combine R1 and R2 alignment hits from paired reads to the 51 functional marker gene prokaryotic database, aggregate to marker gene level and transform to TPM  

# load libraries
library(stringr)
library(tidyr)

# first get sample specific / total counts / total rpk, calculated based on KEGG alignments
setwd(as.character(args[1]))
all_sample_info_files <- list.files()
all_sample_info <- list()

for(f in 1:length(all_sample_info_files)){
  all_sample_info[[f]] <- read.csv(all_sample_info_files[f], sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F)
}
names(all_sample_info) <- gsub(".txt", "", all_sample_info_files)

# read in the the DIAMOND MetMarkDB alignment files and process
setwd(args[2])
all_files <- list.files()

# just select the R1 files
all_files <- all_files[grep("_R1_", all_files)]

# save the results here
agg_hits_list <- list()

# do in a loop where one sample is processed at each iteration of the loop
keep_items <- ls()
keep_items <- c(keep_items, "keep_items")
for(f in 1:length(all_files)){
  setwd(args[2])
  
  file_r1 <- all_files[f]
  file_r2 <- gsub("_R1_", "_R2_", all_files[f])
  
  # the file name format is fixed and defined by the shell script running the DIAMOND alignments
  sample_namn <- gsub("_R1_metmarkdb.txt", "", file_r1)
  print(paste("Processing sample:", sample_namn))
  
  # read-in the input files, this format is fixed and follows the one defined in the shell script for the DIAMOND search run.
  # DIAMOND results for R1
  r1_func <- read.csv(file_r1, sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F, 
                      colClasses = c("character","character","numeric", "integer", "integer", "integer", 
                                     "integer", "integer","numeric","numeric", "numeric", "numeric", "integer"))
  colnames(r1_func) <- c("qseqid", "stitle", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp", "slen")
  print("Read 1 file:")
  print(file_r1)
  
  # DIAMOND results for R2
  r2_func <- read.csv(file_r2, sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F, 
                      colClasses = c("character","character","numeric", "integer", "integer", "integer", 
                                     "integer", "integer","numeric","numeric", "numeric", "numeric", "integer"))
  colnames(r2_func) <- c("qseqid", "stitle", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp", "slen")
  print("Read 2 file:")
  print(file_r2)
  
  # what was imported?
  print(paste("Read-in R1 file with", nrow(r1_func), "hits"))
  print(paste("Read-in R2 file with", nrow(r2_func), "hits"))
  
  # an extra option to filter the alignment hits with more stringent criteria
  if(as.numeric(args[4]) == 1){
    print("Filtering hits with evalue <= 10^-6 and bitscore >= 50")
    r1_func <- r1_func[which(r1_func$evalue <= 10^(-6) & r1_func$bitscore >= 50),]
    r2_func <- r2_func[which(r2_func$evalue <= 10^(-6) & r2_func$bitscore >= 50),]
    
    print(paste("R1 hits remaining after filtering:", nrow(r1_func)))
    print(paste("R2 hits remaining after filtering:", nrow(r2_func)))
  }
  
  # parse input data
  print("Parsing input data and doing RPK transformations")
  
  # parse gene level information.
  r1_genes <- lapply(r1_func[,2], function(x) strsplit(x = as.character(x[1]), split = "-"))
  r1_genes <- unlist(lapply(r1_genes, function(x) x[[1]][1]))
  
  r2_genes <- lapply(r2_func[,2], function(x) strsplit(x = as.character(x[1]), split = "-"))
  r2_genes <- unlist(lapply(r2_genes, function(x) x[[1]][1]))
  
  # join gene information with other
  r1_func <- cbind(r1_genes, r1_func)
  r2_func <- cbind(r2_genes, r2_func)
  
  # for each read pair, we want to calculate hits to exactly the same gene (or protein) only once. So if both read pairs have a hit in the same gene, remove pair2 hit.
  
  # transform to RPK already at this stage. 
  r1_counts <- rep(1, nrow(r1_func))
  r1_func <- cbind(r1_func, r1_counts)
  
  # for length, use amino acid sequence length multiplied by 3 (as 3 nucleotides make a codon) and divide by 1000 (to get kilobases)
  r1_rpk <- as.numeric(r1_func$r1_counts)/((as.numeric(r1_func$slen)*3)/1000)
  r1_func <- cbind(r1_func, r1_rpk)
  
  # R2
  r2_counts <- rep(1, nrow(r2_func))
  r2_func <- cbind(r2_func, r2_counts)
  
  r2_rpk <- as.numeric(r2_func$r2_counts)/((as.numeric(r2_func$slen)*3)/1000)
  r2_func <- cbind(r2_func, r2_rpk)
  
  colnames(r1_func)[1:2] <- c("gene", "read_location")
  colnames(r2_func)[1:2] <- c("gene", "read_location")
  
  # first filter both read pair hits separately according to limits suggested by the creators of the database for short reads
  # then combine the remaining hits for the paired reads.
  
  # filter the marker genes with the limits given by Pok Leung. The DIAMOND aligment should be performed also with the
  # same parameters (now applied) he has given in his scripts for this filtering to be effective.
  # from him: I usually only analyse reads from 130 - 150 bps and retain hits with over 80% query coverage. Rare metabolic marker genes should still be captured when the number of reads is over 5M. 
  # his scripts as well as the functional 51 marker genes protein database has been downloaded from google drive
  # https://drive.google.com/drive/folders/1IeDF4iUbTtPhuaBOyQkrAlXO67EbkJ-A
  
  # a marine microbiome study using similar (but not exactly the same) scripts (but same filtering limits)
  # https://doi.org/10.1038/s41564-023-01322-0

  # reference for the database
  # Leung, Pok Man; Greening, Chris (2021). Compiled Greening lab metabolic marker gene databases. Monash University. Online resource. https://doi.org/10.26180/14431208.v1
  # https://doi.org/10.26180/14431208.v1
  
  print("Filtering marker gene hits according to given threshold for percentage of identical matches")
  # functional genes, categories and limits extracted from the scripts by Pok Leung
  marker_genes <- c("DsrA-", "FCC-", "Sqr-", "Sor-", "AsrA-", "SoxB-", "RbcL-", "AcsB-", "CooS-", "AclB-", "Mcr-", "HbsC-", "HbsT-",
                    "AmoA-", "NxrA-", "NarG-", "NapA-", "NirS-", "NirK-", "NrfA-", "NosZ-", "HzsA-", "NifH-", "NorB-", "Nod-",
                    "CoxL-", "FeFe-", "[Fe]", "NiFe-", "PmoA-", "MmoA-", "McrA-", "IsoA-", "RHO-", "PsaA-", "PsbA-", "ArsC-",
                    "MtrB-","OmcB-","RdhA-", "YgfK-", "ARO-", "Cyc2-", "FdhA-", "SdhA_FrdA-", "AtpA-", "CcoN-", "CoxA-",
                    "CydA-", "CyoA-", "NuoF-")
  
  # pident limits for the different functional genes copied from their scripts and used in https://doi.org/10.1038/s41564-023-01322-0
  marker_gene_limits <- c(50,50,50,50,50,50,60,50,50,50,50,50,75,60,60,50,50,50,50,50,50,50,50,50,50,60,60,50,50,50,60,50,70,50,80,70,50,50,50,50,70,70,50,50,50,70,50,50,50,50,60)
  
  # marker gene categories similarly
  marker_gene_categories <- c(rep("Sulfur cycle",6), rep("Carbon fixation",7), rep("Nitrogen cycle",12), rep("Trace gas metabolism",8), rep("Phototrophy",3),
                              rep("Alternative_e_acceptor",5), rep("Alternative_e_donor",3), rep("Respiration",7))
  
  # remove all marker gene hits with pident (=percentage of identical matches) less than given limit / threshold
  
  # R1 hits
  print("Filtering R1 hits")
  all_rem_locs <- numeric()
  for(i in 1:length(marker_genes)){
    
    # grep all the hits for the specific genes
    locs_int <- grep(marker_genes[i], r1_func$stitle, fixed = T)
    
    if(length(locs_int)>0){
      temp <- r1_func[locs_int,]
      
      # remove those marker gene hits with pident < threshold
      rem_locs <- which(as.numeric(temp$pident) < marker_gene_limits[i])
      
      if(length(rem_locs)>0){
        
        if(i==29){ # NiFe-, copied from their scripts, further filter Group4 hits differently
          temp2 <- temp[-rem_locs,]
          group4_locs <- grep("Group 4", temp2$stitle, fixed = T)
          
          rem_locs1 <- locs_int[rem_locs]
          rem_locs2 <- which(as.numeric(temp2[group4_locs,4])<60)
          
          if(length(rem_locs2)>0){
            locs_int <- locs_int[-rem_locs]
            rem_locs2 <- locs_int[group4_locs[rem_locs2]]
            
            all_rem_locs=c(all_rem_locs, rem_locs1, rem_locs2)
          }else{
            all_rem_locs=c(all_rem_locs, rem_locs1)
          }
          
          rm(temp2)
          rm(group4_locs)
          rm(rem_locs1)
          rm(rem_locs2)
          
        }else{
          rem_locs <- locs_int[rem_locs]
          all_rem_locs <- c(all_rem_locs, rem_locs)
        }
        
        rm(rem_locs)
        rm(temp)
      }
      rm(locs_int)
    }
  }
  
  # save only those hits above the threshold
  if(length(all_rem_locs)>0){
    r1_func_filtered <- r1_func[-all_rem_locs,]
  }else{
    r1_func_filtered <- r1_func
  }
  
  print(paste("Succesfully filtered", length(all_rem_locs), "hits not fulfilling pident criteria from R1 hits" ))
  print(paste(nrow(r1_func_filtered), "hits from a total of", nrow(r1_func), "R1 hits remaining"))
  rm(all_rem_locs)
  
  # R2 hits
  print("Filtering R2 hits")
  all_rem_locs <- numeric()
  for(i in 1:length(marker_genes)){
    
    # grep all the hits for the specific genes
    locs_int <- grep(marker_genes[i], r2_func$stitle, fixed = T)
    
    if(length(locs_int)>0){
      temp <- r2_func[locs_int,]
      
      # remove those marker gene hits with pident < threshold
      rem_locs <- which(as.numeric(temp$pident)<marker_gene_limits[i])
      
      if(length(rem_locs)>0){
        
        if(i==29){ # NiFe-, copiend from their scripts, further filter Group4 hits differently
          temp2 <- temp[-rem_locs,]
          group4_locs <- grep("Group 4", temp2$stitle, fixed = T)
          
          rem_locs1 <- locs_int[rem_locs]
          
          rem_locs2 <- which(as.numeric(temp2[group4_locs,4])<60)
          if(length(rem_locs2)>0){
            locs_int <- locs_int[-rem_locs]
            rem_locs2 <- locs_int[group4_locs[rem_locs2]]
            
            all_rem_locs <- c(all_rem_locs, rem_locs1, rem_locs2)
          }else{
            all_rem_locs <- c(all_rem_locs, rem_locs1)
          }
          
          rm(temp2)
          rm(group4_locs)
          rm(rem_locs1)
          rm(rem_locs2)
          
        }else{
          rem_locs <- locs_int[rem_locs]
          all_rem_locs <- c(all_rem_locs, rem_locs)
        }
        
        rm(rem_locs)
        rm(temp)
      }
      rm(locs_int)
    }
  }
  
  if(length(all_rem_locs)>0){
    r2_func_filtered <- r2_func[-all_rem_locs,]
  }else{
    r2_func_filtered <- r2_func
  }
  
  print(paste("Succesfully filtered", length(all_rem_locs), "hits not fulfilling pident criteria from R2 hits" ))
  print(paste(nrow(r2_func_filtered), "hits from a total of", nrow(r2_func), "R2 hits remaining"))
  rm(all_rem_locs)
  
  # combine hits from the paired reads
  print("Combining the hits from the paired reads")
  int_reads <- intersect(r1_func_filtered$read_location, r2_func_filtered$read_location)
  rownames(r1_func_filtered) <- r1_func_filtered$read_location
  rownames(r2_func_filtered) <- r2_func_filtered$read_location
  
  # read pairs for which hits have been found in both read pairs of the same read
  r1_func_int <- r1_func_filtered[int_reads,]
  r2_func_int <- r2_func_filtered[int_reads,]
  
  # unique read pairs, for which there are only hits in R2 - always keep these
  r2_uniq <- r2_func_filtered[-which(r2_func_filtered$read_location%in%int_reads),]
  
  # order based on read pair location for R1 and R2 - so in the same order for both reads
  r1_func_int <- r1_func_int[order(r1_func_int$read_location),]
  r2_func_int <- r2_func_int[order(r2_func_int$read_location),]
  
  # this works, because there is only one (best) hit for each read and they are in the same order in the mate pairs
  # so works only with one hit / read in each read pair. If more hits / read are included this scripts needs to be modified
  
  # where the R2 hit has been to the same metabolic marker gene. Do not count these twice.
  locs_rem <- which(r2_func_int$gene==r1_func_int$gene) 
  
  # remove these duplicated gene hits for R2 and keep only those truly different hits
  r2_keep <- r2_func_int[-locs_rem,]
  
  if(nrow(r2_keep)>0){
    r2_func_filt <- rbind(r2_uniq, r2_keep)
  }else{
    r2_func_filt <- r2_uniq
  }
  # here we have the R2 hits with a different gene or not hit for the R1 read. So then all gene hits, regardless of did it originate from R1, R2 or both, are counted only once.
  
  print(paste("Preserved all of R1 hits with:", nrow(r1_func_filtered), "hits"))
  print("Number of R2 hits that were not found with paired reads in R1:")
  print(nrow(r2_func_filt))
  print("Making a proportion of all R2 filtered hits:")
  print(nrow(r2_func_filt)/nrow(r2_func_filtered))
  
  # add information for each hit from which read it came from
  r1_func_filtered <- cbind(r1_func_filtered, rep("R1",nrow(r1_func_filtered)))
  r2_func_filt <- cbind(r2_func_filt, rep("R2", nrow(r2_func_filt)))
  
  colnames(r1_func_filtered)[ncol(r1_func_filtered)] <- "hit_readpair"
  colnames(r2_func_filt)[ncol(r2_func_filt)] <- "hit_readpair"
  
  # combine all into one
  colnames(r1_func_filtered) <- colnames(r2_func_filt)
  all_hits <- rbind(r1_func_filtered, r2_func_filt)
  colnames(all_hits)[15:16] <- c("count", "rpk")
  
  print(paste("Succesfully combined R1 and R2 with", nrow(all_hits), "hits in total"))
  
  # summarize the hits into sequence level
  print("Aggregating to sequence level")
  all_hits_aggregated <- aggregate(cbind(count, rpk) ~ stitle, all_hits, sum)
  colnames(all_hits_aggregated)=c("sequence", "count", "rpk")
  all_hits_aggregated$gene=sub("-.*", "", all_hits_aggregated$sequence)
  all_hits_aggregated$sample <- sample_namn
  print(paste("Succesfully aggregated to", nrow(all_hits_aggregated), "unique sequences" ))
  
  # save the combined hits
  setwd(args[3])
  print("Saving the processed files")
  
  dir.create("combined_hits")
  setwd("combined_hits")
  fil_n <- paste(sample_namn, ".txt", sep = "")
  write.table(x = all_hits, file = fil_n, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  setwd("../")
  
  # save the filtered hits for R1 and R2
  dir.create("filtered_hits")
  setwd("filtered_hits")
  fil_n <- paste(sample_namn, "_R1.txt", sep = "")
  write.table(x = r1_func_filtered, file = fil_n, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  fil_n <- paste(sample_namn, "_R2.txt", sep = "")
  write.table(x = r2_func_filt, file = fil_n, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  setwd("../")
  
  # save the combined hits aggregated to sequnce level
  dir.create("aggregated_sequence_hits")
  setwd("aggregated_sequence_hits")
  fil_n <- paste(sample_namn, ".txt", sep = "")
  write.table(x = all_hits_aggregated, file = fil_n, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  setwd("../")
  
  # save the combined hits aggregated to the sequence level also in a list for downstream processing
  agg_hits_list[[f]] <- all_hits_aggregated
  print("Succesfully saved everything")
  
  print("****************************************")
  print("****************************************")
  print(paste(round(f/length(all_files), 3), "of total files processed"))
  print("****************************************")
  print("****************************************")
  
  del_items <- ls()
  del_items <- del_items[-which(del_items%in%keep_items)]
  rm(list = del_items)
}

# name the list
names(agg_hits_list) <- gsub("_R1_metmarkdb.txt","",all_files)

print("Succesfully saved")
print("****************************************")
print("****************************************")

# combine all the hits from the different samples
print("Combining hits from the diffeent samples")
hits_all_samples <- do.call("rbind", agg_hits_list)
print("All samples combined")
print("Summarizing the detected hits to gene level and converting into a matrix form for RPK and counts. Calculating TPM from the RPK.")

# summarize RPK and counts by gene x sample
hits_summary_rpk <- aggregate(rpk ~ sample+gene, hits_all_samples, sum)
hits_summary_count <- aggregate(count ~ sample+gene, hits_all_samples, sum)

# put into a gene x sample matrix
hits_samples_rpk <- spread(hits_summary_rpk, sample, rpk)
rownames(hits_samples_rpk) <- hits_samples_rpk$gene
hits_samples_rpk <- hits_samples_rpk[,-which(colnames(hits_samples_rpk)%in%"gene")]

hits_samples_count <- spread(hits_summary_count, sample, count)
rownames(hits_samples_count) <- hits_samples_count$gene
hits_samples_count <- hits_samples_count[,-which(colnames(hits_samples_count)%in%"gene")]

# set NA to zeroes, change this if needed
hits_samples_rpk[is.na(hits_samples_rpk)] <- 0
hits_samples_count[is.na(hits_samples_count)] <- 0

# now, transform the RPK to TPM by using the RPK scaling factors for the samples from the KEGG alignment. 
# the KEGG alignment is a fairly good estimate of all the prokaryotic genes in the samples and thus the total gene counts and total RPK in the samples.

# similarly, transform to reads per million, RPM, data by using the sample specific count per million scaling factors
rpk_scaling_factors <- numeric(length(all_sample_info))
count_scaling_factors <- numeric(length(all_sample_info))
for(i in 1:length(all_sample_info)){
  temp <- all_sample_info[[i]]
  rownames(temp) <- as.character(temp[,1])
  rpk_scaling_factors[i] <- as.numeric(temp["RPK_PerMillion_Scaling_Factor",2])
  count_scaling_factors[i] <- as.numeric(temp["Count_PerMillion_Scaling_Factor",2])
}
names(rpk_scaling_factors) <- names(all_sample_info)
names(count_scaling_factors) <- names(all_sample_info)

# make sure the order of the samples is correct - similar to the aggregated data
rpk_scaling_factors <- rpk_scaling_factors[colnames(hits_samples_rpk)]
count_scaling_factors <- count_scaling_factors[colnames(hits_samples_count)]

# to get TPM, divide the RPK for each marker gene with the RPK per million scaling factor for each sample
hits_samples_tpm <- hits_samples_rpk
for(i in 1:ncol(hits_samples_tpm)){hits_samples_tpm[,i] <- hits_samples_tpm[,i]/rpk_scaling_factors[i]}

# get RPM similarly
hits_samples_rpm <- hits_samples_count
for(i in 1:ncol(hits_samples_rpm)){hits_samples_rpm[,i] <- hits_samples_rpm[,i]/count_scaling_factors[i]}

# these matrices can be taken forward for analysis.
print("Saving final matrices")

setwd(args[3])
write.table(x = hits_all_samples, file = "Combined_Filtered_Hits_All_Samples.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = hits_summary_rpk, file = "RPK_by_Samples_Genes.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(x = hits_samples_rpk, file = "RPK_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(x = hits_samples_rpm, file = "RPM_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(x = hits_samples_tpm, file = "TPM_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(x = hits_summary_count, file = "Counts_by_Samples_Genes.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(x = hits_samples_count, file = "Counts_Samplewise_Matrix.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

save(hits_samples_rpk, hits_samples_rpm, hits_samples_tpm, hits_samples_count, all_sample_info, rpk_scaling_factors, count_scaling_factors, file = "All_The_Output_Matrices_For_Downstream.RData")

# print out session info
print("SessionInfo:")
sessionInfo()

