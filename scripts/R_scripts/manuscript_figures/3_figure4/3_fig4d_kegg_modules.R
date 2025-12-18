# a small script to plot the abundance and expression of the selectec core KEGG metabolic modules

# load libraries
library(pheatmap)
library(RColorBrewer)

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load datas and put into a list
datas_list <- list()

# metagenomics
load(paste(project_root, "/metagenomics/kegg_diamond/Matrices_For_Downstream.RData", sep = ""))
datas_list[[1]] <-  ko_tpm_data

# metatranscriptomics
load(paste(project_root, "/metatranscriptomics/kegg_diamond/Matrices_For_Downstream.RData", sep = ""))
datas_list[[2]] <-  ko_tpm_data
names(datas_list) <- c("kegg_mg", "kegg_mt")

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# retrieve and parse the downloaded KEGG annotation for KO groups
ko_annotations <- read.csv(file = paste(project_root, "/metadata/2025_02_28_KO_Name_Mapping_REST.txt", sep = "")
                           , header = F, fill = F, sep = "\t", stringsAsFactors = F, quote = "")

all_kos <- as.character(ko_annotations[,1])
ko_annotations <- data.frame(ko_annotations[,-1])
rownames(ko_annotations) <- all_kos
colnames(ko_annotations) <- c("name")

# only include the selected interesting modules for the enrichment analysis
modules_curated <- read.csv(file = paste(project_root, "/metadata/selected_KEGG_modules.txt", sep = ""), quote = "", fill=F, stringsAsFactors = F, sep = "\t", header = F)

# get module identifiers
ko_modules_curated  <- strsplit(x = modules_curated[,1], split = " ")
ko_modules_curated <- unlist(lapply(ko_modules_curated, function(x) x[1]))

# read-in the module completeness information for samples from anvi'o
# set working directory - needs to exist
setwd(paste(project_root, "/downstream/kegg_modules", sep = ""))

# list to save the module information
present_modules <- list()

keep_items <- ls()
keep_items <- c(keep_items, "keep_items")
for(dat in 1:length(datas_list)){
  
  # save module presence in data.frame
  module_presence <- data.frame(matrix(nrow = length(ko_modules_curated), ncol = ncol(datas_list[[dat]])))
  rownames(module_presence) <- ko_modules_curated
  
  setwd(names(datas_list)[dat])
  
  # get the list of samples
  all_samples <- list.files()
  
  for(sam in 1:length(all_samples)){
    setwd(all_samples[sam])
    
    # read-in the module completeness estimates from anvi'o
    modules_info <- read.csv("metabolism_modules.txt", sep = "\t")
    
    # select only the modules complete enough - 75% currently
    modules_3_4 <- modules_info[which(modules_info$stepwise_module_is_complete=="True"),]
    
    # include only the complete enough modules from the modules of interest - curated modules
    modules_selected <- modules_3_4[which(modules_3_4$module%in%ko_modules_curated),]
    
    # save presence information
    module_presence[modules_selected$module,sam] <- 1
    
    setwd("../")
  }
  
  # save sample names
  colnames(module_presence) <- all_samples
  
  # remove those modules which are missing in all samples
  module_presence <- module_presence[rowSums(is.na(module_presence))!=ncol(module_presence), ]
  
  # save
  present_modules[[dat]] <- module_presence
  
  setwd("../")
  
  # clean temporary variables created inside the loop, keep datas_list / modules / present_modules
  del_items <- ls()
  del_items <- del_items[-which(del_items%in%keep_items)]
  rm(list = del_items)
}
names(present_modules) <- names(datas_list)

# separate into mg and mt
modules_data_mg <- present_modules[[1]]
modules_data_mt <- present_modules[[2]]

# get the ko groups in the modules
setwd(paste(project_root, "/downstream/kegg_modules/kegg_mg/P1", sep = ""))

modules_info <- read.csv("metabolism_modules.txt", sep = "\t")
rownames(modules_info) <- modules_info$module
modules_info_mg <- modules_info[rownames(modules_data_mg),]
modules_info_mt <- modules_info[rownames(modules_data_mt),]

# metagenomics
ko_modules_mg <- list()
for(j in 1:nrow(modules_info_mg)){
  ko_modules_mg[[j]] <- unique(unlist(strsplit(x = gsub("[-]", ",",gsub("[+]", ",",gsub(" ", ",",gsub("[)]", "", gsub("[(]","",modules_info_mg$module_definition[j]))))), split = ",")))
}
names(ko_modules_mg) <- paste(modules_info_mg$module, modules_info_mg$module_name, sep = "_")

# metatranscriptomics
ko_modules_mt <- list()
for(j in 1:nrow(modules_info_mt)){
  ko_modules_mt[[j]] <- unique(unlist(strsplit(x = gsub("[-]", ",",gsub("[+]", ",",gsub(" ", ",",gsub("[)]", "", gsub("[(]","",modules_info_mt$module_definition[j]))))), split = ",")))
}
names(ko_modules_mt) <- paste(modules_info_mt$module, modules_info_mt$module_name, sep = "_")

## process metagenomics data
# load the rpk data and the library sizes for samples
load(paste(project_root, "/metagenomics/kegg_diamond/Matrices_For_Downstream.RData", sep = ""))

rpk_mg <- ko_samples_rpk
colnames(rpk_mg) <- paste("P",unlist(lapply(colnames(rpk_mg), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
rpk_mg <- rpk_mg[,rownames(metadata)]
rpk_scaling_factors <- rpk_scaling_factors[rownames(metadata)]

# use only those ko groups included in the low expression filtered tpm data data
rpk_mg <- rpk_mg[rownames(ko_tpm_data),]

# make module data and correct the library sizes with the sum rpk of multicounted ko groups for each sample
multi_rpk <- numeric(ncol(rpk_mg))

# prepare the data frame, where to store the summarized rpk for the modules
module_heat_mg <- modules_data_mg

# organize according to the same sample order
module_heat_mg <- module_heat_mg[,colnames(rpk_mg)]

# go through the samples
for(j in 1:ncol(module_heat_mg)){
  
  # save which ko goups are used
  locs_int_all <- numeric()
  
  # go through the modules
  for(i in 1:nrow(module_heat_mg)){
    # if the module is not present in the sample, skip - NAs remain NAs
    if(is.na(module_heat_mg[i,j])){next}
    
    # otherwise summarize the rpk of all the ko groups in that module for that sample
    locs_int <- which(rownames(rpk_mg)%in%ko_modules_mg[[i]])
    module_heat_mg[i,j] <- sum(rpk_mg[locs_int,j], na.rm = T)
    locs_int_all <- c(locs_int_all, locs_int)
  }
  
  # get multicounted ko groups
  if(any(duplicated(locs_int_all))){
    
    # create a mapping for the multicounts
    
    # get the multicounted kos
    all_ko_mg <- rownames(rpk_mg)[locs_int_all]
    
    # unique kos in this
    unique_ko_mg <- unique(all_ko_mg)
    
    # count how many this each ko is mapped
    ko_multimapping_mg <- unlist(lapply(unique_ko_mg, function(x) length(which(all_ko_mg%in%x))))
    names(ko_multimapping_mg) <- unique_ko_mg
    
    # remove the first count - not a multicount
    ko_multimapping_mg <- ko_multimapping_mg-1
    
    # remove those KO groups not multicounted
    ko_multimapping_mg <- ko_multimapping_mg[-which(ko_multimapping_mg==0)]
    
    # here we have the KO groups that were counted more than once and how many times they were counted more than once
    # correct the library size estimated with this information - may vary accross samples - retain sample comparability
    # although most likely not that many and does not create a big difference
    multi_rpk[j] <- sum(rpk_mg[names(ko_multimapping_mg),j]*ko_multimapping_mg)
  }
}
names(multi_rpk) <- colnames(rpk_mg)

# correct the rpk scaling factors to account for the multicounting

# turn the original rpk scaling factors into samplewise rpk sums
rpk_scaling_factors_multi <- rpk_scaling_factors * 10 ^6

# add the multi rpk for each sample
rpk_scaling_factors_multi <- rpk_scaling_factors_multi + multi_rpk

# make into tpm scaling factors again
rpk_scaling_factors_multi <- rpk_scaling_factors_multi / 10 ^ 6 

# make into tpm data
module_heat_mg_tpm <- module_heat_mg
for(i in 1:ncol(module_heat_mg_tpm)){module_heat_mg_tpm[,i] <- module_heat_mg_tpm[,i]/rpk_scaling_factors_multi[i]}

# save the scaling factors in a variable
rpk_scaling_factors_multi_mg <- rpk_scaling_factors_multi

# save also the multicount rpk in a variable
multi_rpk_mg <- multi_rpk

## process metatranscriptomics data similarly
# load the rpk data and the library sizes for samples
load(paste(project_root, "/metatranscriptomics/kegg_diamond/Matrices_For_Downstream.RData", sep = ""))
rpk_mt <- ko_samples_rpk

# parse column names
colnames(rpk_mt) <- paste("P",unlist(lapply(colnames(rpk_mt), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

# organize into the correct order according metadata
rpk_mt <- rpk_mt[,rownames(metadata)]

# make sure the rpk scaling factors are in the same order
rpk_scaling_factors <- rpk_scaling_factors[rownames(metadata)]

# use only those ko groups included in the low expression filtered TPM data data
rpk_mt <- rpk_mt[rownames(ko_tpm_data),]

# make module data and correct the library sizes with the sum rpk of multicounted ko groups for each sample
multi_rpk <- numeric(ncol(rpk_mt))

# prepare the data frame, where to store the summarized rpk for the modules
module_heat_mt <- modules_data_mt

# organize according to the same sample order
module_heat_mt <- module_heat_mt[,colnames(rpk_mt)]

# go through the samples
for(j in 1:ncol(module_heat_mt)){
  
  # save which ko goups are used
  locs_int_all <- numeric()
  
  # go through the modules
  for(i in 1:nrow(module_heat_mt)){
    # if the module is not present in the sample, skip - NAs remain NAs
    if(is.na(module_heat_mt[i,j])){next}
    
    # otherwise summarize the rpk of all the ko groups in that module for that sample
    locs_int <- which(rownames(rpk_mt)%in%ko_modules_mt[[i]])
    module_heat_mt[i,j] <- sum(rpk_mt[locs_int,j], na.rm = T)
    locs_int_all <- c(locs_int_all, locs_int)
  }
  
  # get multicounted ko groups
  if(any(duplicated(locs_int_all))){
    
    # create a mapping for the multicounts
    
    # get the multicounted kos
    all_ko_mt <- rownames(rpk_mt)[locs_int_all]
    
    # unique kos in this
    unique_ko_mt <- unique(all_ko_mt)
    
    # count how many this each ko is mapped
    ko_multimapping_mt <- unlist(lapply(unique_ko_mt, function(x) length(which(all_ko_mt%in%x))))
    names(ko_multimapping_mt) <- unique_ko_mt
    
    # remove the first count - not a multicount
    ko_multimapping_mt <- ko_multimapping_mt-1
    
    # remove those ko groups not multicounted
    ko_multimapping_mt <- ko_multimapping_mt[-which(ko_multimapping_mt==0)]
    
    # here we have the ko groups that were counted more than once and how many times they were counted more than once
    # correct the library size estimated with this information - may vary accross samples - retain sample comparability
    # although most likely not that many and does not create a big difference
    multi_rpk[j] <- sum(rpk_mt[names(ko_multimapping_mt),j]*ko_multimapping_mt)
  }
}
names(multi_rpk) <- colnames(rpk_mt)

# correct the rpk scaling factors to account for the multicounting

# turn the original rpk scaling factors into samplewise rpk sums
rpk_scaling_factors_multi <- rpk_scaling_factors * 10 ^6

# add the multi rpk for each sample
rpk_scaling_factors_multi <- rpk_scaling_factors_multi + multi_rpk

# make into tpm scaling factors again
rpk_scaling_factors_multi <- rpk_scaling_factors_multi / 10 ^ 6 

# make into tpm data
module_heat_mt_tpm <- module_heat_mt
for(i in 1:ncol(module_heat_mt_tpm)){module_heat_mt_tpm[,i] <- module_heat_mt_tpm[,i]/rpk_scaling_factors_multi[i]}

# save the scaling factors in a variable
rpk_scaling_factors_multi_mt <- rpk_scaling_factors_multi

# save also the multicount rpk in a variable
multi_rpk_mt <- multi_rpk

# here we have the kegg module values we can use for plotting
# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure4", sep = ""))

# annotations for plotting - this now data or study specific - modify accordingly if needed
annotation_col <- meta_use
for(j in 1:ncol(annotation_col)){annotation_col[,j] <- as.character(annotation_col[,j])}

annotation_colors <- list(
  grazing = c("grazed"="tan1", "ungrazed"="forestgreen"),
  veg_clusters = c("t_ces"="darkolivegreen2", "c_ros"="purple", "c_cho"="deeppink"),
  treatment = c("-S"="dodgerblue2", "CTL"="gainsboro", "+S"="firebrick3")
)

# plot both omics similarly, without any clustering

# all modules in MT present in MG, but not all modules present in MG are present in MT
# thus plot all modules present in MG

# plotting data frame for MG
data_heat_mg <- module_heat_mg_tpm

# prepare a similar data.frame for MT
data_heat_mt <- data.frame(matrix(nrow = nrow(data_heat_mg), ncol = ncol(data_heat_mg)))
rownames(data_heat_mt) <- rownames(data_heat_mg)
colnames(data_heat_mt) <- colnames(data_heat_mg)
data_heat_mt[rownames(module_heat_mt_tpm),] <- module_heat_mt_tpm

# add module information to rownames
rownames(data_heat_mg) <- paste(rownames(data_heat_mg), modules_info[rownames(data_heat_mg),3])
rownames(data_heat_mt) <- paste(rownames(data_heat_mt), modules_info[rownames(data_heat_mt),3])

# arrange according to MG expression, from highest to lowest
row_m_mg <- rowMeans(data_heat_mg, na.rm = T)
row_m_mg <- row_m_mg[order(row_m_mg, decreasing = T)]

# order both omics data.frames similarly
data_heat_mg <- data_heat_mg[names(row_m_mg),]
data_heat_mt <- data_heat_mt[names(row_m_mg),]

# plot in the same plot
r <- as.numeric(unlist(data_heat_mg))
r <- c(quantile(r, 0.05, na.rm = T), quantile(r, 0.95, na.rm = T))
breakslist_mg <- seq(from=r[1], to=r[2], length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))

d1 <- data_heat_mg
d2 <- data_heat_mt

# add omics information to column names
colnames(d1) <- paste(colnames(d1), "MG")
colnames(d2) <- paste(colnames(d2), "MT")

# create common data.frame
common_data <- cbind(d1,d2)

# replicate annotation - same samples
a1 <- annotation_col
a2 <- annotation_col

# add omics information to rownames of annotation frames - will match common data column names
rownames(a1) <- paste(rownames(a1), "MG")
rownames(a2) <- paste(rownames(a2), "MT")

# combine into annotation
common_annotation <- rbind(a1, a2)

# finally plot
pheatmap(mat = common_data, na_col = "darkgrey", fontsize = 12, cluster_rows = F, gaps_col = 36, breaks = breakslist_mg, color = cols_heat, cluster_cols = F, annotation_col = common_annotation, annotation_colors = annotation_colors, scale="none", cellwidth = 14, cellheight = 14, silent = T, fontsize_row = 14, fontsize_col = 14, filename = "KEGG_Modules_both_read.pdf")

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures
