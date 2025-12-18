# a small script to save the metadata regarding vegetation clusters in the study

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load metadata
load(paste(project_root, "/metadata/Study_Metadata.RData", sep = ""))

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

# save as table
meta_table <- rbind(table(meta_use$grazing, meta_use$veg_clusters), table(meta_use$treatment, meta_use$veg_clusters))
write.csv(x = meta_table, file = "metadata_vegetation.csv", quote = F, row.names = T, col.names = T, sep = "\t")

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in excel