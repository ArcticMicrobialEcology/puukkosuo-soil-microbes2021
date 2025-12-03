args <- commandArgs(trailingOnly = TRUE)
# the first input parameter needs to be the metagenomics directory for the phyloflash result csv files
# the second input parameter needs to be the metatranscriptomics directory for the phyloflash result csv files
# the third input parameter needs to be the file for the sample information and study metadata (e.g. Oulanka_ACAP_study_site.xlsx)

# a small script compose taxonomic count data from phyloflash results

# process metagenomics and metatranscriptomics data simultaneously with the same script
for(a in seq(1:2)){
  keep_files <- ls()
  setwd(args[a])
  
  # get the phyloflas output csv files
  phylo_input_files <- list.files()
  
  # use the full abundance (full taxonomic rank) files
  phylo_input_files <- phylo_input_files[grep("full_abundance", phylo_input_files)]
  
  # parse 
  for(i in 1:length(phylo_input_files)){
    fileName <- phylo_input_files[i]
    fileData <- read.csv(fileName,header=F)
    colnames(fileData) <- c("NTU", phylo_input_files[i])
    
    if (!exists("NTUcounts")) {
      NTUcounts <- fileData
    } else {
      NTUcounts <- merge(x=NTUcounts, y=fileData, by="NTU", all=TRUE)
    }
  }
  
  ntu_names      <- NTUcounts$NTU
  sample_names   <- gsub(".phyloFlash.NTUfull_abundance.csv", "", colnames(NTUcounts[,-1]))
  
  NTUcounts  <- NTUcounts[,-1]
  rownames(NTUcounts) <- ntu_names
  colnames(NTUcounts) <- sample_names
  NTUcounts[is.na(NTUcounts)] <- 0     # turn NA into 0
  
  # get taxonomic table
  # summarize to different taxonomical levels / ranks
  taxons <- strsplit(x = rownames(NTUcounts), split = ";")
  
  # for those taxons identified to a certain level only, replace the lower levels with unknown
  # set unknown taxon annotation as unknown instead of the previous higher level as by phyloflash default
  taxons_parsed <- lapply(taxons, function(x) {
    if(any(grepl("[(]", substring(x, 1, 1)))){x[which(grepl("[(]", substring(x, 1, 1)))]="unknown"} 
    x
  })
  
  # turn into a dataframe
  taxons_parsed_df <- data.frame(t(as.data.frame(taxons_parsed)), check.names = F, stringsAsFactors = F)
  colnames(taxons_parsed_df) <- c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(taxons_parsed_df) <- rownames(NTUcounts)
  # this is not correct for Eukaryotes - but don't care about them - at least for now
  
  # set all taxa containing keywords such as metagenome, uncultured into unknown for reliability
  key_w <- c("metagenome", "uncultured")
  
  for(i in 1:ncol(taxons_parsed_df)){
    for(j in 1:length(key_w)){
      if(any(grepl(key_w[j], taxons_parsed_df[,i]))){
        locs_int <- grep(key_w[j], taxons_parsed_df[,i])
        taxons_parsed_df[locs_int, i] <- "unknown"
      }
    }
  }
  
  # load study metadata
  load(as.character(args[3]))

  # parse tables
  otu_table_all <- NTUcounts
  tax_table_all <-taxons_parsed_df

  # summarize Eukaryotic counts
  euk_locs <-  which(tax_table_all$Domain%in%"Eukaryota") # these can be just summarized as one
  euk_counts <- colSums(otu_table_all[euk_locs,])
  
  # make otu and tax tables with just one line for Eukaryotes
  otu_table <- otu_table_all[-euk_locs,]
  tax_table <- tax_table_all[-euk_locs,]
  
  otu_table <- rbind(otu_table, euk_counts)
  tax_table <- rbind(tax_table, c("Eukaryota", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown"))
  rownames(otu_table)[nrow(otu_table)] <- rownames(tax_table)[nrow(tax_table)] <- "Eukaryota;(Eukaryota);(Eukaryota);(Eukaryota);(Eukaryota);(Eukaryota);(Eukaryota)"
  
  # parse sample names to match metadata
  colnames(otu_table) <- paste("P",unlist(lapply(colnames(otu_table), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
  colnames(otu_table_all) <- paste("P",unlist(lapply(colnames(otu_table_all), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
  
  # organize according to metadata
  rownames(metadata) <- metadata$`Short code`
  otu_table <- otu_table[,rownames(metadata)]
  otu_table_all <- otu_table_all[,rownames(metadata)]
  
  # save everything
  setwd("../")
  dir.create("downstream")
  setwd("downstream")
  
  # make blocks unique, block 1 in grazed is not related to block 1 in ungrazed.
  metadata$Block[which(metadata$Grazing=="ungrazed")] <- as.numeric(metadata$Block[which(metadata$Grazing=="ungrazed")])+6
  
  # convert sample_table variables to factors
  for(i in 1:ncol(metadata)){
    metadata[,i] <- as.factor(metadata[,i])
  }
  
  save(otu_table, otu_table_all, tax_table, tax_table_all, metadata, file = "Otu_Tax_Tables_Parsed.RData")
  
  del_files <- ls()
  del_files <- del_files[-which(del_files%in%keep_files)]
  rm(list = del_files)
}

# print out session info
print("SessionInfo:")
sessionInfo()