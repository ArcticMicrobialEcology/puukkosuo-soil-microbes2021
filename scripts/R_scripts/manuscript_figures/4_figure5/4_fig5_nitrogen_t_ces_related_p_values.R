# a small script to prepare some environmental variables for use for the later plotting and analysis scripts.

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load parsed gas flux and pore water data
load(paste(project_root, "/metadata/environmental_data/Meth_fluxes_porewater.RData", sep = ""))

# load metadata
load(paste(project_root, "/metadata/Study_Metadata.RData", sep = ""))

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# separate grazing into two classes based on vegetation
meta_use$grazing_vegetation <- factor(paste(meta_use$grazing, meta_use$veg_clusters, sep = "_"), levels=c("ungrazed_t_ces","grazed_t_ces", "ungrazed_c_cho", "grazed_c_cho", "grazed_c_ros"))

# pore water variables
# loook at pore water variables in the grazing / exclusion treatment in only the T.ces vegetation cluster
# nitrate + nitrite
temp_data <- cbind(log(pore_water_nitrogen[,1]), meta_use)
colnames(temp_data)[1] <- "nitrate_nitrite"

# compare only T.ces plots grazed vs. ungrazed
temp_data_tces <- temp_data[which(temp_data$veg_clusters=="t_ces"),]

# wilcoxon test
p_val_t_ces_nitrate_wilcox <- wilcox.test(x = as.numeric(temp_data_tces$nitrate_nitrite[which(temp_data_tces$grazing=="grazed")]), y=as.numeric(temp_data_tces$nitrate_nitrite[which(temp_data_tces$grazing=="ungrazed")]),
                                          alternative = "two.sided")
# p-value = 0.0318

# total nitrogen
temp_data <- cbind(log(pore_water_nitrogen[,2]), meta_use)
colnames(temp_data)[1] <- "total_nit"

# compare only T.ces plots grazed vs. ungrazed
temp_data_tces <- temp_data[which(temp_data$veg_clusters=="t_ces"),]

# wilcoxon test
p_val_t_ces_total_nit_wilcox <- wilcox.test(x = as.numeric(temp_data_tces$total_nit[which(temp_data_tces$grazing=="grazed")]), y=as.numeric(temp_data_tces$total_nit[which(temp_data_tces$grazing=="ungrazed")]),
                                            alternative = "two.sided")
# p-value = 0.09133

# lookt at also gene ratios in MG and MT data
# MG
# load data
load(paste(project_root, "/metagenomics/metmarkdb_diamond/Matrices_For_Downstream.RData", sep = ""))
data <- metmark_tpm_data

# nirK+nirS/nosZ)
nit_ratio2 <- (as.numeric(data["NirK",]) +  as.numeric(data["NirS",])) / as.numeric(data["NosZ",])
temp_data <- cbind(nit_ratio2, meta_use)
colnames(temp_data) <- c("nit_ratio2", colnames(meta_use))
temp_data_tces <- temp_data[which(temp_data$veg_clusters=="t_ces"),]

# wilcoxon test
p_val_t_ces_nit_ratio2_mg_wilcox <- wilcox.test(x = as.numeric(temp_data_tces$nit_ratio2[which(temp_data_tces$grazing=="grazed")]), y=as.numeric(temp_data_tces$nit_ratio2[which(temp_data_tces$grazing=="ungrazed")]),
                                                alternative = "two.sided")
# p-value = 0.4442

# MT
load(paste(project_root, "/metatranscriptomics/metmarkdb_diamond/Matrices_For_Downstream.RData", sep = ""))
data <- metmark_tpm_data

# nirK+nirS/nosZ)
nit_ratio2 <- (as.numeric(data["NirK",]) +  as.numeric(data["NirS",])) / as.numeric(data["NosZ",])
temp_data <- cbind(nit_ratio2, meta_use)
colnames(temp_data) <- c("nit_ratio2", colnames(meta_use))
temp_data_tces <- temp_data[which(temp_data$veg_clusters=="t_ces"),]

# wilcoxon test
p_val_t_ces_nit_ratio2_mt_wilcox <- wilcox.test(x = as.numeric(temp_data_tces$nit_ratio2[which(temp_data_tces$grazing=="grazed")]), y=as.numeric(temp_data_tces$nit_ratio2[which(temp_data_tces$grazing=="ungrazed")]),
                                                alternative = "two.sided")
# p-value = 0.06218

# save p-values
ses_info <- sessionInfo()

# save into figure5 directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure5", sep = ""))

save(p_val_t_ces_nitrate_wilcox, p_val_t_ces_total_nit_wilcox,
     p_val_t_ces_nit_ratio2_mg_wilcox, p_val_t_ces_nit_ratio2_mt_wilcox, ses_info, file = "Nitrogen_T_Ces_P_values.RData")

# print out session info
print("SessionInfo:")
sessionInfo()

