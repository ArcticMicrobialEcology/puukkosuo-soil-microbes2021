# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load the datas and metadata
load(paste(project_root, "/metadata/Study_Metadata.RData", sep = ""))

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))
rownames(meta_use) <- metadata$`Short code`

# get condition and sample names, initialize directories
all_conditions <- as.character(read.table(paste(project_root, "/metadata/Chosen_Conditions.txt", sep = ""))[,1]) 
sample_info_dir <- paste(project_root, "/metadata", sep = "")
all_sample_names <- as.character(read.table(paste(project_root, "/metadata/MG_Sample_Names.txt", sep = ""))[,1]) 

# get bowtie mapping stats
bowtie_dir <- paste(project_root, "/metagenomics/megahit_assemblies/bowtie_evaluation", sep = "")

# gather bowtie overall alignment rates for samples
map_stats <- data.frame(matrix(nrow = length(all_sample_names), ncol = length(all_conditions)))
rownames(map_stats) <- all_sample_names
colnames(map_stats) <- c(all_conditions)

for(i in 1:length(all_conditions)){
  
  # get sample names for condition
  setwd(sample_info_dir)
  sample_names <- as.character(read.table(paste("MG_Sample_Names_", all_conditions[i],".txt", sep = ""))[,1])
  setwd(bowtie_dir)
  
  # gather stats
  cond_stats <- numeric(length(sample_names))
  for(j in 1:length(sample_names)){
    
    # read-in the bowtie log file
    file_n <- paste(all_conditions[i],"_",sample_names[j],"_bowtie_log.txt", sep = "")
    temp_file <- readLines(file_n)
    
    # save the overall alignment rate
    locs_int <- grep("overall alignment rate", temp_file)
    al_stats <- temp_file[locs_int]
    al_stats <- as.numeric(gsub("% overall alignment rate", "", al_stats))
    cond_stats[j] <- al_stats
    
    rm(file_n)
    rm(temp_file)
    rm(locs_int)
    rm(al_stats)
  }
  names(cond_stats) <- sample_names
  map_stats[names(cond_stats),i] <- cond_stats
  
  rm(sample_names)
  rm(cond_stats)
}

# parse sample names
rownames(map_stats) <- paste("P",unlist(lapply(rownames(map_stats), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

# arrange according to metadata
map_stats <- map_stats[rownames(meta_use),]

# put into one numerical variable
map_stats_numeric <- numeric(nrow(map_stats))
map_stats_numeric[which(meta_use$grazing=="grazed")] <- as.numeric(map_stats$Grazed[which(meta_use$grazing=="grazed")])
map_stats_numeric[which(meta_use$grazing=="ungrazed")] <- as.numeric(map_stats$UnGrazed[which(meta_use$grazing=="ungrazed")])
names(map_stats_numeric) <- rownames(map_stats)

# plot the gathered statistics
# change directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

cols <- c(rep("tan1", length(which(metadata$Grazing=="grazed"))), rep("forestgreen", length(which(metadata$Grazing=="ungrazed"))))
{
  pdf(file = "coassemblies_overall_alignment_rate.pdf", width = 21, height = 7)
  par(mar = c(10.1, 4.1, 4.1, 2.1))
  barplot(map_stats_numeric, beside = T,las=2, col=cols, ylab = "Overall alignment rate %")
  legend("topleft", legend = c("Outside", "Inside"), fill = c("tan1", "forestgreen"))
  dev.off()
}

# save as tables also
write.csv(x = map_stats, file = "coassemblies_overall_alignment_rate.csv")

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures
