# a small script to plot the differentially abundant (MG) and differentially expressed KO groups (MT) from the lmm analysis

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(pheatmap)
library(RColorBrewer)
library(openxlsx)

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
annot_use <- ko_annotations

# load the lmm analysis results
# metagenomics
load(paste(project_root, "/downstream/lmm_analysis/kegg_mg/lmm_analysis_done.RData", sep = ""))

# pick the relevant data
data <- datas_list$kegg_mg

# define the DE lists
# grazing
no_interaction_features <- rownames(sig_values_lmm_best)[-which(sig_values_lmm_best$best_model=="interaction")]
grazing_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(1)], sig_values_lmm_best[no_interaction_features,c(1)], fdr_strict_best[no_interaction_features,c(1)]) 
colnames(grazing_res_simple) <- c("coef", "pvalue", "fdr")
rownames(grazing_res_simple) <- no_interaction_features
de_grazing_simple <- grazing_res_simple[which(grazing_res_simple$fdr<=0.1),]
de_grazing_simple <- de_grazing_simple[order(de_grazing_simple$pvalue),]
temp <- de_grazing_simple
# 19 DA ko groups

# add annotations
annots_plot <- character(nrow(temp))
for(ro in 1:nrow(temp)){annots_plot[ro] <-  paste(c(rownames(temp)[ro], as.character(annot_use[rownames(temp)[ro],])), collapse = "; ")}
rownames(temp) <- annots_plot

# save table into excel
# create a workbook
wb <- createWorkbook()

# add worksheet
addWorksheet(wb, "DA_grazing_MG")
writeData(wb, sheet = "DA_grazing_MG", x = temp, startRow = 1, startCol = 1, rowNames = TRUE)

# snow increase
snow_increase_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(2)], sig_values_lmm_best[no_interaction_features,c(2)], fdr_strict_best[no_interaction_features,c(2)]) 
colnames(snow_increase_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_increase_res_simple) <- no_interaction_features
de_snow_increase_simple <- snow_increase_res_simple[which(snow_increase_res_simple$fdr<=0.1),]
de_snow_increase_simple <- de_snow_increase_simple[order(de_snow_increase_simple$pvalue),]
# nothing detected

# snow decrease
snow_decrease_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(3)], sig_values_lmm_best[no_interaction_features,c(3)], fdr_strict_best[no_interaction_features,c(3)]) 
colnames(snow_decrease_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_decrease_res_simple) <- no_interaction_features
de_snow_decrease_simple <- snow_decrease_res_simple[which(snow_decrease_res_simple$fdr<=0.1),]
de_snow_decrease_simple <- de_snow_decrease_simple[order(de_snow_decrease_simple$pvalue),]
# 1 ko group

temp <- de_snow_decrease_simple

# add annotations
annots_plot <- character(nrow(temp))
for(ro in 1:nrow(temp)){annots_plot[ro] <-  paste(c(rownames(temp)[ro], as.character(annot_use[rownames(temp)[ro],])), collapse = "; ")}
rownames(temp) <- annots_plot

# save table into excel
# add worksheet
addWorksheet(wb, "DA_snow_decrease_MG")
writeData(wb, sheet = "DA_snow_decrease_MG", x = temp, startRow = 1, startCol = 1, rowNames = TRUE)

# ko groups with interactions btw. exclusion and snow treatments
interaction_features <- rownames(sig_values_lmm_best)[which(sig_values_lmm_best$best_model=="interaction")]
grazing_treatment_interaction <- data.frame(coefficients_lmm_best[interaction_features,c(4:5)], sig_values_lmm_best[interaction_features,c(4:5)], fdr_strict_best[interaction_features,c(4:5)], stringsAsFactors = F, check.names = F) 
colnames(grazing_treatment_interaction) <- paste(c("coef", "coef", "pvalue", "pvalue", "fdr", "fdr"), colnames(grazing_treatment_interaction))
rownames(grazing_treatment_interaction) <- interaction_features
de_grazing_treatment_interaction <- grazing_treatment_interaction[which(grazing_treatment_interaction$`fdr grazinggrazed:treatment+S`<=0.1 | grazing_treatment_interaction$`fdr grazinggrazed:treatment-S`<=0.1),]
# one da ko group
temp <- de_grazing_treatment_interaction

# add annotations
annots_plot <- character(nrow(temp))
for(ro in 1:nrow(temp)){annots_plot[ro] <-  paste(c(rownames(temp)[ro], as.character(annot_use[rownames(temp)[ro],])), collapse = "; ")}
rownames(temp) <- annots_plot

# save table into excel
# add worksheet
addWorksheet(wb, "DA_grazing_snow_int_MG")
writeData(wb, sheet = "DA_grazing_snow_int_MG", x = temp, startRow = 1, startCol = 1, rowNames = TRUE)

# start plotting
# grazing as heatmap
data_heat <- data[rownames(de_grazing_simple),]

# add annotations
annots_plot <- character(nrow(data_heat))
for(ro in 1:nrow(data_heat)){annots_plot[ro] <-  paste(c(rownames(data_heat)[ro], as.character(annot_use[rownames(data_heat)[ro],])), collapse = "; ")}
rownames(data_heat) <- annots_plot

# define scale and colors
breakslist <- seq(from=-2, to=2, length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))

# annotations for plotting
annotation_col <- meta_use
for(j in 1:ncol(annotation_col)){annotation_col[,j] <- as.character(annotation_col[,j])}

annotation_colors <- list(
  grazing = c("grazed"="tan1", "ungrazed"="forestgreen"),
  veg_clusters = c("t_ces"="darkolivegreen2", "c_ros"="purple", "c_cho"="deeppink"),
  treatment = c("-S"="dodgerblue2", "CTL"="gainsboro", "+S"="firebrick3")
)


# plot
# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure4", sep = ""))
pheatmap(mat = data_heat, cluster_rows = TRUE, cluster_cols = FALSE, scale="row", annotation_col = annotation_col, annotation_colors = annotation_colors, gaps_col = 18, cellwidth = 14, cellheight = 14, breaks=breakslist, color = cols_heat, silent = T, fontsize_row = 8, fontsize_col = 8, filename = "de_grazing_KEGG_MG_read.pdf", main = "")

# snow decrease - only 1 DE feature, plot boxplot
vals <- as.numeric(data[rownames(de_snow_decrease_simple),])

# add annotations
r_name <- rownames(de_snow_decrease_simple)
main_text <-  paste(c(r_name, as.character(annot_use[r_name,])), collapse = "; ")

# wrap
main_text <- paste(strwrap(main_text, width = 20), collapse = "\n")

# plot - add some empty plots for easier postprocessing in inkscape
{
  pdf(file = "de_snow_decrease_KEGG_MG.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plot(1,1)
  plot(1,1)
  g <- as.factor(meta_use$treatment)
  r_boxplot <- range(vals, na.rm = T)
  boxplot(vals~g, col=c("gainsboro", "firebrick3", "dodgerblue2"), xlab = "", ylim = r_boxplot, ylab = "Copies per million", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  dev.off()
}

# grazing snow treatment interaction
# snow decrease - only 1 DE feature, plot boxplot
vals <- as.numeric(data[rownames(de_grazing_treatment_interaction),])

# add annotations
r_name <- rownames(de_grazing_treatment_interaction)
main_text <-  paste(c(r_name, as.character(annot_use[r_name,])), collapse = "; ")

# wrap
main_text <- paste(strwrap(main_text, width = 20), collapse = "\n")

# plot - add some empty plots for easier postprocessing in inkscape
{
  pdf(file = "de_treatment_grazing_KEGG_MG_read.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plot(1,1)
  
  g <- as.factor(paste(meta_use$treatment, meta_use$grazing, sep = " "))
  g <- factor(g, levels = c("CTL ungrazed", "CTL grazed", "+S ungrazed", "+S grazed", "-S ungrazed", "-S grazed"))
  r_boxplot <- range(vals, na.rm = T)
  boxplot(vals~g, col=c("forestgreen", "tan1"), las=2, xlab = "", ylim = r_boxplot, ylab = "Copies per million", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  dev.off()
}

# metatranscriptomics
load(paste(project_root, "/downstream/lmm_analysis/kegg_mt/lmm_analysis_done.RData", sep = ""))

# pick the relevant data
data <- datas_list$kegg_mt

# define the DE lists
# grazing
no_interaction_features <- rownames(sig_values_lmm_best)[-which(sig_values_lmm_best$best_model=="interaction")]
grazing_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(1)], sig_values_lmm_best[no_interaction_features,c(1)], fdr_strict_best[no_interaction_features,c(1)]) 
colnames(grazing_res_simple) <- c("coef", "pvalue", "fdr")
rownames(grazing_res_simple) <- no_interaction_features
de_grazing_simple <- grazing_res_simple[which(grazing_res_simple$fdr<=0.1),]
de_grazing_simple <- de_grazing_simple[order(de_grazing_simple$pvalue),]
temp <- de_grazing_simple
# one de ko group

# add annotations
annots_plot <- character(nrow(temp))
for(ro in 1:nrow(temp)){annots_plot[ro] <-  paste(c(rownames(temp)[ro], as.character(annot_use[rownames(temp)[ro],])), collapse = "; ")}
rownames(temp) <- annots_plot

# save table into excel
# add worksheet
addWorksheet(wb, "DE_grazing_MT")
writeData(wb, sheet = "DE_grazing_MT", x = temp, startRow = 1, startCol = 1, rowNames = TRUE)

# snow increase
snow_increase_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(2)], sig_values_lmm_best[no_interaction_features,c(2)], fdr_strict_best[no_interaction_features,c(2)]) 
colnames(snow_increase_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_increase_res_simple) <- no_interaction_features
de_snow_increase_simple <- snow_increase_res_simple[which(snow_increase_res_simple$fdr<=0.1),]
de_snow_increase_simple <- de_snow_increase_simple[order(de_snow_increase_simple$pvalue),]
# nothing

# snow decrease
snow_decrease_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(3)], sig_values_lmm_best[no_interaction_features,c(3)], fdr_strict_best[no_interaction_features,c(3)]) 
colnames(snow_decrease_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_decrease_res_simple) <- no_interaction_features
de_snow_decrease_simple <- snow_decrease_res_simple[which(snow_decrease_res_simple$fdr<=0.1),]
de_snow_decrease_simple <- de_snow_decrease_simple[order(de_snow_decrease_simple$pvalue),]
temp <- de_snow_decrease_simple
# one de ko group

# add annotations
annots_plot <- character(nrow(temp))
for(ro in 1:nrow(temp)){annots_plot[ro] <-  paste(c(rownames(temp)[ro], as.character(annot_use[rownames(temp)[ro],])), collapse = "; ")}
rownames(temp) <- annots_plot

# save table into excel
# add worksheet
addWorksheet(wb, "DE_snow_decrease_MT")
writeData(wb, sheet = "DE_snow_decrease_MT", x = temp, startRow = 1, startCol = 1, rowNames = TRUE)

# ko groups with interactions btw. exclusion and snow treatments
interaction_features <- rownames(sig_values_lmm_best)[which(sig_values_lmm_best$best_model=="interaction")]
grazing_treatment_interaction <- data.frame(coefficients_lmm_best[interaction_features,c(4:5)], sig_values_lmm_best[interaction_features,c(4:5)], fdr_strict_best[interaction_features,c(4:5)], stringsAsFactors = F, check.names = F) 
colnames(grazing_treatment_interaction) <- paste(c("coef", "coef", "pvalue", "pvalue", "fdr", "fdr"), colnames(grazing_treatment_interaction))
rownames(grazing_treatment_interaction) <- interaction_features
de_grazing_treatment_interaction <- grazing_treatment_interaction[which(grazing_treatment_interaction$`fdr grazinggrazed:treatment+S`<=0.1 | grazing_treatment_interaction$`fdr grazinggrazed:treatment-S`<=0.1),]
# nothing

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure4", sep = ""))
# save compiled excel workbook
saveWorkbook(wb, "figure4_KO_DA_DE_taxa.xlsx", overwrite = TRUE)

# plot
# grazing - only 1 DE feature, plot boxplot
vals <- as.numeric(data[rownames(de_grazing_simple),])

# add annotations
r_name <- rownames(de_grazing_simple)
main_text <-  paste(c(r_name, as.character(annot_use[r_name,])), collapse = "; ")

# wrap
main_text <- paste(strwrap(main_text, width = 20), collapse = "\n")

# plot - add empty plots for easier postprocessing in inkscape
{
  pdf(file = "de_grazing_KEGG_MT.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plot(1,1)
  plot(1,1)
  g <- as.factor(meta_use$grazing)
  r_boxplot <- range(vals, na.rm = T)
  boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, ylab = "Copies per million", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  dev.off()
}

# snow decrease - only 1 DE feature, plot boxplot
vals <- as.numeric(data[rownames(de_snow_decrease_simple),])

# add annotations
r_name <- rownames(de_snow_decrease_simple)
main_text <-  paste(c(r_name, as.character(annot_use[r_name,])), collapse = "; ")

# wrap
main_text <- paste(strwrap(main_text, width = 20), collapse = "\n")

# plot
{
  pdf(file = "de_snow_decrease_KEGG_MT.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plot(1,1)
  plot(1,1)
  g <- as.factor(meta_use$treatment)
  r_boxplot <- range(vals, na.rm = T)
  boxplot(vals~g, col=c("gainsboro", "firebrick3", "dodgerblue2"), xlab = "", ylim = r_boxplot, ylab = "Copies per million", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  dev.off()
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures

