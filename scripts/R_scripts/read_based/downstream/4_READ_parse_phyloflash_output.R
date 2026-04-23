args <- commandArgs(trailingOnly = TRUE)
# the first input parameter needs to be the metagenomics directory for the phyloflash result csv files
# the second input parameter needs to be the metatranscriptomics directory for the phyloflash result csv files
# the third input parameter needs to be the RData file for the sample information and study metadata (e.g. information in Oulanka_ACAP_study_site.xlsx)
# the fourth input parameter needs to be the directory for the silva database mapping files

# a small script compose taxonomic count data from phyloflash results

# read in the silva ssu taxonomy mapping files, the file names need to be in the standard format given when downloading the files from silva
setwd(as.character(args[4]))

# taxonomy map with organims names
silva_ssu_taxmap <- read.table("taxmap_slv_ssu_ref_138.1.txt", sep = "\t", fill=T, quote = "", comment.char = "", header = T)

# taxonomy rank information
silva_ssu_ranks <- read.table("tax_slv_ssu_138.1.txt", sep = "\t", fill=T, quote = "", comment.char = "", header = F)
silva_ssu_ranks <- silva_ssu_ranks[,-4] # empty column
colnames(silva_ssu_ranks) <- c("path", "taxid", "rank", "version")

# process metagenomics and metatranscriptomics data simultaneously with the same script
for(a in seq(1:2)){

  print(paste("Processing data from the directory:", as.character(args[a])))

  keep_files <- ls()
  setwd(as.character(args[a]))
  
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
  
  # try to gather the proper ranks for each hit into a table
  # get all possible ranks in the silva ssu database
  all_unique_ranks <- unique(silva_ssu_ranks$rank)
  taxons_parsed_df <- data.frame(matrix(nrow = length(taxons_parsed), ncol = length(all_unique_ranks)+1))
  colnames(taxons_parsed_df) <- c(all_unique_ranks, "organism_name")
  
  # a variable to to check consistency if custom / funny organism names are included in the output
  # e.g. Bacteria;Elusimicrobiota;Lineage IIc;Elusimicrobia bacterium RIFOXYA2_FULL_39_19;
  na_checks <- logical(length(taxons_parsed))
  
  # start mapping
  for(i in 1:length(taxons_parsed)){
    temp_taxon <- taxons_parsed[[i]]
    
    # parse the whole string, remove unknown
    if(any(temp_taxon=="unknown")){
      temp_taxon_stripped <- temp_taxon[-which(temp_taxon=="unknown")]
    }else{
      temp_taxon_stripped <- temp_taxon
    }
    
    # try to get ranks
    ranks_taxon <- character(length(temp_taxon_stripped))
    for(j in 1:length(temp_taxon_stripped)){
      temp_path <- paste(paste(temp_taxon_stripped[1:j], collapse = ";"), ";", sep = "")
      
      # do we find this directly in the silva ssu rank file
      if(any(silva_ssu_ranks$path%in%temp_path)){
        loc_rank <- which(silva_ssu_ranks$path%in%temp_path)
        ranks_taxon[j] <- silva_ssu_ranks$rank[loc_rank]
      }else{
        ranks_taxon[j] <- NA
      }
    }
    
    # if the last rank is na, e.g. it is a custom / funny name, see if we can find this in the mapping file
    if(is.na(ranks_taxon[length(ranks_taxon)])){
      
      # the ranks to be used, discarding the last entry, the custom / funny name
      temp_taxon_save <- temp_taxon_stripped[-length(temp_taxon_stripped)]
      names(temp_taxon_save) <- ranks_taxon[-length(ranks_taxon)]
      
      # get the organism name as the last entry
      temp_name <- temp_taxon_stripped[length(temp_taxon_stripped)]
      
      # save it just in case
      taxons_parsed_df$organism_name[i] <- temp_name
      
      # make a new path without the last entry
      temp_path <- paste(paste(temp_taxon_stripped[1:(length(temp_taxon_stripped)-1)], collapse = ";"), ";", sep = "")
      
      # do we find this name in the silva ssu taxmap file
      if(any(silva_ssu_taxmap$organism_name%in%temp_name)){
        
        # there might be several
        loc_name <- which(silva_ssu_taxmap$organism_name%in%temp_name)
        
        # subset with the name
        temp_silva <- silva_ssu_taxmap[loc_name,]
        
        # see if any paths match to the organism path examined
        path_check <- any(temp_silva$path==temp_path)
        
        if(path_check){
          
          # get the taxid for the matching paths
          temp_taxid <- temp_silva$taxid[which(temp_silva$path==temp_path)]
          
          if(length(temp_taxid)>0){
            
            # match this to the rank file
            loc_rank_taxid <- which(silva_ssu_ranks$taxid%in%temp_taxid)
            
            # get only the unique ranks, these should all be the same
            temp_rank_taxid <- unique(silva_ssu_ranks$rank[loc_rank_taxid])
            
            # check that this matches the last rank in the obtained ranks
            if(ranks_taxon[length(ranks_taxon)-1]==temp_rank_taxid){
              na_checks[i] <- TRUE
            }
          } else{
            na_checks[i] <- FALSE
          }
        } else{
          na_checks[i] <- FALSE
        }
        
      }else{
        na_checks[i] <- FALSE
      }
    } else{
      na_checks[i] <- NA
      temp_taxon_save <- temp_taxon_stripped
      names(temp_taxon_save) <- ranks_taxon
    }
    
    # for some eukaryotes there might be multiple annotations at the same rank level, e.g. major_clade
    # in this case, collapse these into one. If processed downstream at those rank levels, must be handled somehow
    if(any(duplicated(names(temp_taxon_save)))){ 
      uniq_taxon_ranks <- unique(names(temp_taxon_save))
      new_taxon_save <- character(length(uniq_taxon_ranks))
      for(r in 1:length(uniq_taxon_ranks)){
        new_taxon_save[r] <- paste(temp_taxon_save[which(names(temp_taxon_save)%in%uniq_taxon_ranks[r])], collapse = ";")
      }
      temp_taxon_save <- new_taxon_save
      names(temp_taxon_save) <- uniq_taxon_ranks
    }
    # save ranks into the table
    taxons_parsed_df[i, c(names(temp_taxon_save))] <- temp_taxon_save
  }
  
  # see if any of the taxons with custom / funny organism names last obtained rank didn't match the rank given in the silva ssu taxonomy rank file
  if(!any(!na_checks, na.rm = T)){
    print("all organism names found in the taxon mapping file")
  } else {
    print("some organism names not found in the taxon mapping file. Please check these.")
  }
  
  # remove all na columns
  na_cols <- apply(taxons_parsed_df, 2, function(x) all(is.na(x)))
  taxons_parsed_df <- taxons_parsed_df[, -c(which(na_cols))]
  rownames(taxons_parsed_df) <- rownames(NTUcounts)
  
  # now since the focus on here is on prokaryotes, summarize all the eukaryotic taxonomic annotations to simplify the table
  taxons_prok <- taxons_parsed_df
  euk_locs <-  which(taxons_prok$domain%in%"Eukaryota") 
  
  euk_frame <- data.frame(matrix(nrow = length(euk_locs), ncol = ncol(taxons_prok)), stringsAsFactors = F, check.names = F)
  colnames(euk_frame) <- colnames(taxons_prok)
  euk_frame$domain <- rep("Eukaryota", nrow(euk_frame))
  rownames(euk_frame) <- rownames(taxons_prok)[euk_locs]
  taxons_prok[euk_locs,] <- euk_frame
  
  # remove all na columns again
  na_cols <- apply(taxons_prok, 2, function(x) all(is.na(x)))
  taxons_prok <- taxons_prok[, -c(which(na_cols))]
  
  # organize in phylogenetically descending order
  
  # main ranks included for prokayotes (bacteria and archeae in silva)
  main_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species", "organism_name")
  main_ranks_included <- main_ranks[which(main_ranks%in%colnames(taxons_prok))]
  taxons_prok <- taxons_prok[,c(main_ranks_included)]
  
  # capitalize first letter for better compatability with downstream scripts
  colnames(taxons_prok) <- paste0(toupper(substring(colnames(taxons_prok), 1, 1)),substring(colnames(taxons_prok), 2))
  
  # replace NAs with unknown
  taxons_prok[is.na(taxons_prok)] <- "unknown"
  
  # set all taxa containing keywords such as metagenome, uncultured into unknown for reliability
  key_w <- c("metagenome", "uncultured")
  
  for(i in 1:ncol(taxons_prok)){
    for(j in 1:length(key_w)){
      if(any(grepl(key_w[j], taxons_prok[,i]))){
        locs_int <- grep(key_w[j], taxons_prok[,i])
        taxons_prok[locs_int, i] <- "unknown"
      }
    }
  }
  
  # do some sanity checks
  print("Are there any taxons with a Genus annotations but not a Domain annotation")
  print(with(taxons_prok, table(Domain == "unknown" & Genus != "unknown")))
  
  print("How many have genus but missing Family?")
  print(table(taxons_prok$Genus != "unknown" & taxons_prok$Family == "unknown"))
  
  print("How many have genus but missing Order?")
  print(table(taxons_prok$Family != "unknown" & taxons_prok$Order  == "unknown"))
  
  print("Do any ranks contain semicolons from the “duplicate-rank collapse”?")
  print(sapply(taxons_prok, function(col) sum(grepl(";", col))))
  
  # load study metadata
  load(as.character(args[3]))
  
  # parse tables
  otu_table_all <- NTUcounts
  tax_table_all <- taxons_prok
  
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
  
  save(otu_table, otu_table_all, tax_table, tax_table_all, taxons_prok, taxons_parsed_df, metadata, file = "Otu_Tax_Tables_Parsed.RData")
  
  del_files <- ls()
  del_files <- del_files[-which(del_files%in%keep_files)]
  rm(list = del_files)
}

# print out session info
print("Session Info:")
sessionInfo()
