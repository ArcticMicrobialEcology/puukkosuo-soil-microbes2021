args <- commandArgs(trailingOnly = TRUE)
# first input parameter needs to be the metadata directory for the study. The metadata directory
# needs to contain the essential information files, such as the Oulanka_ACAP_study_site.xlsx, Puukkosuo_vegetatation_clusters.xlsx

# the second input parameter needs to be the metagenomics raw data directory - from which the sample names are extracted
# the third input parameter needs to be the metatranscriptomics raw data directory - from which the sample names are extracted

# e.g. 
# args <- character(3)
# args[1] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/metadata"
# args[2] <- "/scratch/project_2009164/2_OULANKA/1_RAW/metagenomics"
# args[3] <- "/scratch/project_2009164/2_OULANKA/1_RAW/metatranscriptomics"

# a setup script to setup sample names in the experimental conditions etc. for the analysis. 
# modify as needed

# load libraries
library(readxl)

# read in the metadata
setwd(as.character(args[1]))
metadata <- data.frame(read_excel(path = "Oulanka_ACAP_study_site.xlsx"), stringsAsFactors = F, check.names = F)
rownames(metadata) <- metadata$Plot

# add vegetations clusters, vegetation cluster information received from Eeva Järvi-Laturi.Based on cluster and indicator species analysis.
veg_cluster <- data.frame(read_excel(path = "Puukkosuo_vegetation clusters.xlsx"), stringsAsFactors = F, check.names = F)
rownames(veg_cluster) <- veg_cluster$plot
veg_cluster <- veg_cluster[rownames(metadata),]
veg_cluster$cluster <- gsub("C.Ros", "c_ros", veg_cluster$cluster)
veg_cluster$cluster <- gsub("C.Cho", "c_cho", veg_cluster$cluster)
veg_cluster$cluster <- gsub("T.Ces", "t_ces", veg_cluster$cluster)
metadata$Veg_clusters <- veg_cluster$cluster

# define raw data directories
mg_raw_data_dir <- as.character(args[2])
mt_raw_data_dir <- as.character(args[3])

# get file names
setwd(mg_raw_data_dir)
mg_files <- list.files()

setwd(mt_raw_data_dir)
mt_files <- list.files()

# shorten into unique sample names
mg_samples <- unique(paste(unlist(lapply(mg_files, function(x) strsplit(x = x, split = "L004")[[1]][1])), "L004", sep = ""))
mt_samples <- unique(paste(unlist(lapply(mt_files, function(x) strsplit(x = x, split = "L004")[[1]][1])), "L004", sep = ""))

# save
setwd(as.character(args[1]))
write.table(x = mg_samples, file = "MG_Sample_Names.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples, file = "MT_Sample_Names.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# get sample names related to experimental conditions
# parse sample names in data to match sample names in metadata
names(mg_samples) <- paste("P",unlist(lapply(mg_samples, function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
names(mt_samples) <- paste("P",unlist(lapply(mt_samples, function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

# order according to metadata¨
mg_samples <- mg_samples[metadata$`Short code`]
mt_samples <- mt_samples[metadata$`Short code`]

# exclusion treatment (grazing)
# samples related to the grazing condition
mg_samples_grazing <- mg_samples[which(metadata$Grazing=="grazed")]
mt_samples_grazing <- mt_samples[which(metadata$Grazing=="grazed")]

# samples related to the ungrazing condition
mg_samples_ungrazing <- mg_samples[which(metadata$Grazing=="ungrazed")]
mt_samples_ungrazing <- mt_samples[which(metadata$Grazing=="ungrazed")]

# snow treatment
# samples related to the snow addition treatment
mg_samples_snow_addition <- mg_samples[which(metadata$Treatment=="+S")]
mt_samples_snow_addition <- mt_samples[which(metadata$Treatment=="+S")]

# samples related to the snow removal treatment
mg_samples_snow_removal <- mg_samples[which(metadata$Treatment=="-S")]
mt_samples_snow_removal <- mt_samples[which(metadata$Treatment=="-S")]

# samples related to the snow control treatment
mg_samples_snow_control <- mg_samples[which(metadata$Treatment=="CTL")]
mt_samples_snow_control <- mt_samples[which(metadata$Treatment=="CTL")]

# microsite not used

# vegetation cluster
# samples related to the Carex Rostrata - c_ros
mg_samples_vegcluster_cros <- mg_samples[which(metadata$Veg_clusters=="c_ros")]
mt_samples_vegcluster_cros <- mt_samples[which(metadata$Veg_clusters=="c_ros")]

# samples related to the Carex chordorrhiza - c_cho
mg_samples_vegcluster_ccho <- mg_samples[which(metadata$Veg_clusters=="c_cho")]
mt_samples_vegcluster_ccho <- mt_samples[which(metadata$Veg_clusters=="c_cho")]

# samples related to the Trichophorum cespitosum - t_ces
mg_samples_vegcluster_tces <- mg_samples[which(metadata$Veg_clusters=="t_ces")]
mt_samples_vegcluster_tces <- mt_samples[which(metadata$Veg_clusters=="t_ces")]

# save
# grazing
write.table(x = mg_samples_grazing, file = "MG_Sample_Names_Grazed.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples_grazing, file = "MT_Sample_Names_Grazed.txt", quote = F, sep = "\t", row.names = F, col.names = F)

write.table(x = mg_samples_ungrazing, file = "MG_Sample_Names_UnGrazed.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples_ungrazing, file = "MT_Sample_Names_UnGrazed.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# snow treatment
write.table(x = mg_samples_snow_addition, file = "MG_Sample_Names_SnowAddition.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples_snow_addition, file = "MT_Sample_Names_SnowAddition.txt", quote = F, sep = "\t", row.names = F, col.names = F)

write.table(x = mg_samples_snow_removal, file = "MG_Sample_Names_SnowRemoval.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples_snow_removal, file = "MT_Sample_Names_SnowRemoval.txt", quote = F, sep = "\t", row.names = F, col.names = F)

write.table(x = mg_samples_snow_control, file = "MG_Sample_Names_SnowControl.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples_snow_control, file = "MT_Sample_Names_SnowControl.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# vegetation cluster
write.table(x = mg_samples_vegcluster_cros, file = "MG_Sample_Names_CRos.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples_vegcluster_cros, file = "MT_Sample_Names_CRos.txt", quote = F, sep = "\t", row.names = F, col.names = F)

write.table(x = mg_samples_vegcluster_ccho, file = "MG_Sample_Names_CCho.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples_vegcluster_ccho, file = "MT_Sample_Names_CCho.txt", quote = F, sep = "\t", row.names = F, col.names = F)

write.table(x = mg_samples_vegcluster_tces, file = "MG_Sample_Names_TCes.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = mt_samples_vegcluster_tces, file = "MT_Sample_Names_TCes.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# save all unique condition suffixes
all_conditions <- c("Grazed", "UnGrazed", "SnowAddition", "SnowRemoval", "SnowControl", "CRos", "CCho", "TCes")
write.table(x = all_conditions, file = "All_Conditions.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# save how the different conditions make up treatment groups - or the unique groups for each condition - not organized in any specific way
tr_groups <- data.frame(matrix(nrow = 3, ncol = 3), stringsAsFactors = F)
colnames(tr_groups) <- c("Grazing", "SnowTreatment", "Veg_clusters")
tr_groups$Grazing <- c("Grazed", "UnGrazed", NA) 
tr_groups$SnowTreatment <- c("SnowAddition", "SnowRemoval", "SnowControl") 
tr_groups$Veg_clusters <- c("CRos", "CCho", "TCes") 
write.table(x = tr_groups, file = "Treatment_Groups.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# save the chosen conditions for the coassembly
chosen_conditions <- c("Grazed", "UnGrazed")
write.table(x = chosen_conditions, file = "Chosen_Conditions.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# save the modified metadata file to be used later in downstream analysis
write.table(x = metadata, file = "Study_Metadata.txt", quote = F, sep = "\t", row.names = F, col.names = T)
save(metadata, file = "Study_Metadata.RData")

# print out session info
print("SessionInfo:")
sessionInfo()
