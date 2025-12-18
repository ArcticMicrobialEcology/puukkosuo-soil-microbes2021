args <- commandArgs(trailingOnly = TRUE)

# the first parameter needs to be the prodigal gene and transcript count directory for the contigs
# the second parameter needs to be the directory for the mapping file for the separate coassembly (MAGs) genes to the contig genes
# the third parameter needs to be the directory for the final dereplicated MAGs 
# the fourth parameter needs to be the read-based KEGG metagenomics filtered normalized data list R object
# the fifth parameter needs to be the read-based KEGG metatranscriptomics filtered normalized data list R object

# a small script to get the gene and transcript counts for all dereplicated final MAGs

# the gene counts need to be mapped to concatenated contig genes and these should be used instead of the coassembly gene counts?

# load the concantenated results
setwd(args[1])

# genes
load("gene_counts_annotations.RData")
gene_sample_counts_concatenated <- gene_sample_counts
gene_info_concatenated <- gene_info
gene_annotations_concatenated <- gene_annotations

# transcripts
load("transcript_counts_annotations.RData")
transcript_sample_counts_concatenated <- transcript_sample_counts
transcript_info_concatenated <- transcript_info
transcript_annotations_concatenated <- transcript_annotations

# load mapping file to concatenated gene annotations - gene_info
setwd(as.character(args[2]))
load("Gene_Contig_Mapping_Sep_Coassem_Conc_Coassem.RData")

# change directory to the final MAG directory
setwd(as.character(args[3]))

# load the dereplication report and other stuff
load("Final_MAG_Taxonomy_Mark_Gene.RData")

all_bins <- list.files()
all_bins <- all_bins[grep("cluster_", all_bins)]

wd <- getwd()

bin_gene_counts <- list()
bin_transcript_counts <- list()
bin_annotations <- list()

keep_items <- ls()
keep_items <- c(keep_items, "keep_items")
for(bin in 1:length(all_bins)){
  print(paste("Processing bin/MAG:", all_bins[bin]))
  setwd(all_bins[bin])
  
  all_files <- list.files()
  bin_genes <- read.csv(file = all_files[grep("-gene_calls.txt", all_files)], sep = "\t", header = T, quote = "", fill = F, stringsAsFactors = F)
  rownames(bin_genes) <- as.character(bin_genes$gene_callers_id)
  
  # get the appropriate gene info
  # get the condition of the final MAG
  loc_bin_final_mag <- which(dereplication_report$cluster%in%all_bins[bin])
  condition_temp <- dereplication_report$condition[loc_bin_final_mag]
  gene_info <- mapping_files[condition_temp]
  gene_info <- gene_info[[1]]
  
  # we only want prodigal genes and not rrna genes here - at least not now?
  prod_genes <- which(bin_genes$gene_callers_id%in%gene_info$Geneid)
  bin_genes <- bin_genes[prod_genes,]
  
  # map to concatenated genes
  bin_mapping <- gene_info[rownames(bin_genes),]
  bin_genes$concat_gene_id <- bin_mapping$New_Geneid
  rownames(bin_genes) <- bin_genes$concat_gene_id
  
  bin_annot <- gene_annotations_concatenated[rownames(bin_genes),] # these should match
  bin_gene_count <- gene_sample_counts_concatenated[rownames(bin_genes),]
  bin_transcript_count <- transcript_sample_counts_concatenated[rownames(bin_genes),]
  
  l1 <- bin_genes$stop-bin_genes$start
  l2 <- bin_annot$stop-bin_annot$start
  length_check <- all(l1==l2)
  
  print(paste("Gene lengths match for bin", all_bins[bin], length_check))
  
  # save everything
  write.table(x = bin_annot, file = paste(all_bins[bin],"-prodigal_gene_annotations.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(x = bin_gene_count, file = paste(all_bins[bin],"-prodigal_gene_counts.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(x = bin_transcript_count, file = paste(all_bins[bin],"-prodigal_transcript_counts.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  bin_annotations[[bin]] <- bin_annot
  bin_gene_counts[[bin]] <- bin_gene_count
  bin_transcript_counts[[bin]] <- bin_transcript_count
  
  del_items <- ls()
  del_items <- del_items[-which(del_items%in%keep_items)]
  rm(list = del_items)
  gc()
  setwd(wd)
}
names(bin_gene_counts) <- all_bins
names(bin_transcript_counts) <- all_bins
names(bin_annotations) <- all_bins

setwd(wd)
save(bin_annotations, bin_gene_counts, bin_transcript_counts, file = "MAG_Gene_Transcript_Counts_Annotations.RData")

# load read-based filtered datas to get the sample metadata and all sample info
# metagenomics
load(as.character(args[4]))
all_sample_info_mg <- all_sample_info

# metatranscriptomics
load(as.character(args[5]))
all_sample_info_mt <- all_sample_info

# store the composed datas here
bin_mg_datas <- list()
bin_mt_datas <- list()

# compile marker gene level tpm, rpk, etc. datas for all final MAGs
keep_items1 <- ls()
keep_items1 <- c(keep_items1, "keep_items1")
for(bin in 1:length(bin_annotations)){
  
  print(paste("Processing bin/MAG:", names(bin_annotations)[bin]))
  
  gene_annotations <- bin_annotations[[bin]]
  
  # parse the MetMarkDB genes
  metmark_db_genes <- strsplit(x = gene_annotations$metmarkdb_function, split = ";")
  metmark_db_genes <- unlist(lapply(metmark_db_genes, function(x) x[1]))
  metmark_db_genes[is.na(metmark_db_genes)]=""
  gene_annotations$metmarkdb_genes <- metmark_db_genes 
  
  # what sort of annotations we want to process - look at only metabolic marker genes for now
  process_cols <- which(colnames(gene_annotations)%in%c("metmarkdb_genes"))
  
  # process both metagenomics and metatranscriptomics similarly
  keep_items2 <- ls()
  keep_items2 <- c(keep_items2, "keep_items2")
  for(om in 1:2){
    
    # change directory
    setwd(wd)
    
    if(om==1){
      gene_sample_counts <- bin_gene_counts[[bin]]
      all_sample_info <- all_sample_info_mg
      print("Processing metagenomics data")
    }else{
      gene_sample_counts <- bin_transcript_counts[[bin]]
      all_sample_info <- all_sample_info_mt
      print("Processing metatranscriptomics data")
    }
    
    # store the results in lists
    count_datas <- list()
    rpk_datas <- list()
    tpm_datas <- list()
    rpk_sample_scale_factors <- list()
    
    for(j in 1:length(process_cols)){
      
      print(paste("Processing", colnames(gene_annotations)[process_cols[j]], "annotations"))
      
      # get all unique features in the data
      all_feat <- unique(unlist(strsplit(x = gene_annotations[,process_cols[j]], split = "; ")))
      if(any(all_feat%in%"")){all_feat <- all_feat[-which(all_feat%in%"")]}
      
      # summarize into feature level
      # rpk
      samples_rpk <- data.frame(matrix(nrow = length(all_feat), ncol = ncol(gene_sample_counts)))
      rownames(samples_rpk) <- all_feat
      colnames(samples_rpk) <- colnames(gene_sample_counts)
      
      # count
      samples_count <- samples_rpk
      
      # get gene lengths for rpk
      # Anvi’o follows the convention of string indexing and splicing that is identical to the way one does it in Python or C. 
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
      
      all_sample_info <- all_sample_info[colnames(samples_rpk)]
      rpk_sample_factors <- numeric(ncol(samples_rpk))
      for(i in 1:length(all_sample_info)){
        temp <- all_sample_info[[i]]
        rownames(temp) <- as.character(temp[,1])
        
        # get the total RPK for each sample
        tot_rpk <- as.numeric(temp["Total_RPK",2])
        
        # calculate the RPK per million scaling factor
        rpk_sample_factors[i] <- tot_rpk / 10^6
      }
      names(rpk_sample_factors) <- colnames(samples_rpk)
      
      
      # to get tpm, divide the rpk for each gene with the rpk per million scaling factor for each sample
      samples_tpm <- samples_rpk
      for(i in 1:ncol(samples_tpm)){samples_tpm[,i] <- samples_tpm[,i]/rpk_sample_factors[i]}
      
      # save all the datas into lists
      count_datas[[j]] <- samples_count
      rpk_datas[[j]] <- samples_rpk
      tpm_datas[[j]] <- samples_tpm
      rpk_sample_scale_factors[[j]] <- rpk_sample_factors
    }
    
    # add suitable names
    namn <- colnames(gene_annotations)[process_cols]
    namn <- gsub("_genes", "", namn)
    names(tpm_datas) <- names(rpk_datas) <- names(count_datas) <- names(rpk_sample_scale_factors) <- namn
    
    # save
    bin_datas <- list(count_datas, rpk_datas, tpm_datas, rpk_sample_scale_factors)
    names(bin_datas) <- c("count", "rpk", "tpm", "rpk_sample_factors")
    
    setwd(names(bin_annotations)[bin])
    
    if(om==1){
      bin_mg_datas[[bin]] <- bin_datas
      save(bin_datas, file = "MG_Count_RPK_TPM_datas.RData")
    }else{
      bin_mt_datas[[bin]] <- bin_datas
      save(bin_datas, file = "MT_Count_RPK_TPM_datas.RData")
    }
    del_items <- ls()
    del_items <- del_items[-which(del_items%in%keep_items2)]
    rm(list = del_items)
  }
  
  del_items <- ls()
  del_items <- del_items[-which(del_items%in%keep_items1)]
  rm(list = del_items)
  
  print("*******************************************************")
  print("*******************************************************")
}
names(bin_mg_datas) <- names(bin_annotations)
names(bin_mt_datas) <- names(bin_annotations)

# save all the datas
setwd(as.character(args[3]))
save(bin_mg_datas, file = "MG_Summarized_Count_Datas_For_All_Bins.RData")
save(bin_mt_datas, file = "MT_Summarized_Count_Datas_For_All_Bins.RData")
