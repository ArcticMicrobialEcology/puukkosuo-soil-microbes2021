# a small script to plot the differentially abundant (MG) and differentially expressed taxa (MT) from the lmm analysis

# load libraries
library(pheatmap)
library(RColorBrewer)
library(phyloseq)
library(microbiome)
library(openxlsx)

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load datas
# metagenomics
load(paste(project_root, "/metagenomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = ""))

# here we have compositional data already aggregated to the order level
psq_order_mg <- psq_order

# get phylum levels as well
psq_phylum_mg <- psq_phylum

# metatranscriptomics 
load(paste(project_root, "/metatranscriptomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = ""))

psq_order_mt <- psq_order

# get phylum levels as well
psq_phylum_mt <- psq_phylum

# get compositional data sets
data_order_mg <- data.frame(otu_table(psq_order_mg), stringsAsFactors = F, check.names = F)
data_order_mt <- data.frame(otu_table(psq_order_mt), stringsAsFactors = F, check.names = F)
data_phylum_mg <- data.frame(otu_table(psq_phylum_mg), stringsAsFactors = F, check.names = F)
data_phylum_mt <- data.frame(otu_table(psq_phylum_mt), stringsAsFactors = F, check.names = F)

# put into a list for plotting
datas_list <- list(data_order_mg, data_phylum_mg, data_order_mt, data_phylum_mt)
names(datas_list) <- c("order_mg", "phylum_mg", "order_mt", "phylum_mt")

# put also the phyloseq objects into a list
psq_list <- list(psq_order_mg, psq_phylum_mg, psq_order_mt, psq_phylum_mt)
names(psq_list) <- c("order_mg", "phylum_mg", "order_mt", "phylum_mt")

# clr transform
psq_clr_list <- list()
datas_clr <- list()

for(i in 1:length(psq_list)){
  temp_psq <- psq_list[[i]]
  temp_psq <-  microbiome::transform(x = temp_psq, transform = "clr")
  
  psq_clr_list[[i]] <- temp_psq
  datas_clr[[i]] <- data.frame(otu_table(temp_psq), stringsAsFactors = F, check.names = F)
}
names(psq_clr_list) <- names(datas_clr) <- names(datas_list)

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# load the results from the lmm analysis
# metagenomics
# order
load(paste(project_root, "/downstream/lmm_analysis/order_mg/lmm_analysis_done.RData", sep = ""))

# transform compositional data into percent proportions
data_order_mg <- data_order_mg * 100

# define the DE lists
# exclusion (grazing here)
no_interaction_features <- rownames(sig_values_lmm_best)[-which(sig_values_lmm_best$best_model=="interaction")]
grazing_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(1)], sig_values_lmm_best[no_interaction_features,c(1)], fdr_strict_best[no_interaction_features,c(1)]) 
colnames(grazing_res_simple) <- c("coef", "pvalue", "fdr")
rownames(grazing_res_simple) <- no_interaction_features
de_grazing_simple <- grazing_res_simple[which(grazing_res_simple$fdr<=0.1),]
de_grazing_simple <- de_grazing_simple[order(de_grazing_simple$pvalue),]
# no significant results

# snow increase
snow_increase_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(2)], sig_values_lmm_best[no_interaction_features,c(2)], fdr_strict_best[no_interaction_features,c(2)]) 
colnames(snow_increase_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_increase_res_simple) <- no_interaction_features
de_snow_increase_simple <- snow_increase_res_simple[which(snow_increase_res_simple$fdr<=0.1),]
de_snow_increase_simple <- de_snow_increase_simple[order(de_snow_increase_simple$pvalue),]
# no significant results

# snow decrease
snow_decrease_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(3)], sig_values_lmm_best[no_interaction_features,c(3)], fdr_strict_best[no_interaction_features,c(3)]) 
colnames(snow_decrease_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_decrease_res_simple) <- no_interaction_features
de_snow_decrease_simple <- snow_decrease_res_simple[which(snow_decrease_res_simple$fdr<=0.1),]
de_snow_decrease_simple <- de_snow_decrease_simple[order(de_snow_decrease_simple$pvalue),]
# no significant results

# features with significant interactions between the exclusion and snow treatments
interaction_features <- rownames(sig_values_lmm_best)[which(sig_values_lmm_best$best_model=="interaction")]
grazing_treatment_interaction <- data.frame(coefficients_lmm_best[interaction_features,c(4:5)], sig_values_lmm_best[interaction_features,c(4:5)], fdr_strict_best[interaction_features,c(4:5)], stringsAsFactors = F, check.names = F) 
colnames(grazing_treatment_interaction) <- paste(c("coef", "coef", "pvalue", "pvalue", "fdr", "fdr"), colnames(grazing_treatment_interaction))
rownames(grazing_treatment_interaction) <- interaction_features
de_grazing_treatment_interaction <- grazing_treatment_interaction[which(grazing_treatment_interaction$`fdr grazinggrazed:treatment+S`<=0.1 | grazing_treatment_interaction$`fdr grazinggrazed:treatment-S`<=0.1),]
# no significant results
# nothing to plot

# phylum level
# load lmm results
load(paste(project_root, "/downstream/lmm_analysis/phylum_mg/lmm_analysis_done.RData", sep = ""))

# transform compositional data into percent proportions
data_phylum_mg <- data_phylum_mg * 100

# define the DE lists
# exclusion (grazing here)
no_interaction_features <- rownames(sig_values_lmm_best)[-which(sig_values_lmm_best$best_model=="interaction")]
grazing_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(1)], sig_values_lmm_best[no_interaction_features,c(1)], fdr_strict_best[no_interaction_features,c(1)]) 
colnames(grazing_res_simple) <- c("coef", "pvalue", "fdr")
rownames(grazing_res_simple) <- no_interaction_features
de_grazing_simple <- grazing_res_simple[which(grazing_res_simple$fdr<=0.1),]
de_grazing_simple <- de_grazing_simple[order(de_grazing_simple$pvalue),]
# no significant results

# snow increase
snow_increase_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(2)], sig_values_lmm_best[no_interaction_features,c(2)], fdr_strict_best[no_interaction_features,c(2)]) 
colnames(snow_increase_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_increase_res_simple) <- no_interaction_features
de_snow_increase_simple <- snow_increase_res_simple[which(snow_increase_res_simple$fdr<=0.1),]
de_snow_increase_simple <- de_snow_increase_simple[order(de_snow_increase_simple$pvalue),]
# no significant results

# snow decrease
snow_decrease_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(3)], sig_values_lmm_best[no_interaction_features,c(3)], fdr_strict_best[no_interaction_features,c(3)]) 
colnames(snow_decrease_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_decrease_res_simple) <- no_interaction_features
de_snow_decrease_simple <- snow_decrease_res_simple[which(snow_decrease_res_simple$fdr<=0.1),]
de_snow_decrease_simple <- de_snow_decrease_simple[order(de_snow_decrease_simple$pvalue),]
# no significant results

# features with significant interactions between the exclusion and snow treatments
interaction_features <- rownames(sig_values_lmm_best)[which(sig_values_lmm_best$best_model=="interaction")]
grazing_treatment_interaction <- data.frame(coefficients_lmm_best[interaction_features,c(4:5)], sig_values_lmm_best[interaction_features,c(4:5)], fdr_strict_best[interaction_features,c(4:5)], stringsAsFactors = F, check.names = F) 
colnames(grazing_treatment_interaction) <- paste(c("coef", "coef", "pvalue", "pvalue", "fdr", "fdr"), colnames(grazing_treatment_interaction))
rownames(grazing_treatment_interaction) <- interaction_features
de_grazing_treatment_interaction <- grazing_treatment_interaction[which(grazing_treatment_interaction$`fdr grazinggrazed:treatment+S`<=0.1 | grazing_treatment_interaction$`fdr grazinggrazed:treatment-S`<=0.1),]
# no significant results

# nothing to plot also here
# metatranscriptomics
# order
# load lmm results
load(paste(project_root, "/downstream/lmm_analysis/order_mt/lmm_analysis_done.RData", sep = ""))

# transform compositional data into percent proportions
data_order_mt <- data_order_mt * 100

# define the DE lists
# grazing
no_interaction_features <- rownames(sig_values_lmm_best)[-which(sig_values_lmm_best$best_model=="interaction")]
grazing_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(1)], sig_values_lmm_best[no_interaction_features,c(1)], fdr_strict_best[no_interaction_features,c(1)]) 
colnames(grazing_res_simple) <- c("coef", "pvalue", "fdr")
rownames(grazing_res_simple) <- no_interaction_features
de_grazing_simple <- grazing_res_simple[which(grazing_res_simple$fdr<=0.1),]
de_grazing_simple <- de_grazing_simple[order(de_grazing_simple$pvalue),]

# save the results table into excel
# create a workbook
wb <- createWorkbook()

# add worksheet
addWorksheet(wb, "DE_grazing_order_MT")

# write the excel
writeData(wb, sheet = "DE_grazing_order_MT", x = de_grazing_simple, startRow = 1, startCol = 1, rowNames = TRUE)

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure3", sep = ""))

# plot heatmap
data <- datas_clr$order_mt
data_heat <- data[rownames(de_grazing_simple),]

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

# add phylum names to plot

# load silva GTDB mapping
load("Silva_GTDB_phylum_mapping.RData")
rownames(silva_gtdb_map) <- silva_gtdb_map$SILVA

# compose order level data with proper taxonomy, as aggregate_rare by microbiome does not preserve taxonomy tables well
psq_order_mt_tax <- aggregate_taxa(x = psq, level = "Order")
tax_mt <- data.frame(tax_table(psq_order_mt_tax), stringsAsFactors = F, check.names = F)
tax_mt <- tax_mt[rownames(data_heat),]

# GTDB phylum names
all_phyla <- tax_mt$Phylum
all_phyla <- silva_gtdb_map[all_phyla, 2]
rownames(data_heat) <- paste(rownames(data_heat), all_phyla, sep = "; ")

# plot heatmap

# italicize names
ph <- pheatmap(mat = data_heat, cluster_rows = TRUE, cluster_cols = FALSE, scale="row", annotation_col = annotation_col, annotation_colors = annotation_colors, gaps_col = 18, cellwidth = 14, cellheight = 16, breaks=breakslist, color = cols_heat, silent = T, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "")

# extract gtable
gt <- ph$gtable

# identify rowname grobs and set fontface = italic
row_ids <- grep("row_names", gt$layout$name)
for (i in row_ids) {
  gr <- gt$grobs[[i]]
  if (inherits(gr, "text")) {
    gr$gp$fontface <- "italic"
    gt$grobs[[i]] <- gr
  }
}

# Export final PDF
pdf("grazing_order_MT_clr_transformed.pdf", width = 12, height = 3)
grid.newpage()
grid.draw(gt)
dev.off()

# snow increase
snow_increase_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(2)], sig_values_lmm_best[no_interaction_features,c(2)], fdr_strict_best[no_interaction_features,c(2)]) 
colnames(snow_increase_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_increase_res_simple) <- no_interaction_features
de_snow_increase_simple <- snow_increase_res_simple[which(snow_increase_res_simple$fdr<=0.1),]
de_snow_increase_simple <- de_snow_increase_simple[order(de_snow_increase_simple$pvalue),]

# save results into excel
# add worksheet
addWorksheet(wb, "DE_snow_increase_order_MT")
writeData(wb, sheet = "DE_snow_increase_order_MT", x = de_snow_increase_simple, startRow = 1, startCol = 1, rowNames = TRUE)

# only one feature - boxplot - use compositional data for the boxplot
vals <- as.numeric(data_order_mt[rownames(de_snow_increase_simple),])

# add annotations
tax_mt <- data.frame(tax_table(psq_order_mt_tax), stringsAsFactors = F, check.names = F)
tax_mt <- tax_mt[rownames(de_snow_increase_simple),]

# GTDB annotations
all_phyla <- tax_mt$Phylum
all_phyla <- silva_gtdb_map[all_phyla, 2]

# main text for the boxplot
main_text <-  paste(rownames(de_snow_increase_simple), all_phyla, sep = "; ")

# wrap
main_text <- paste(strwrap(main_text, width = 20), collapse = "\n")

# italize
main_text <- bquote(italic(.(main_text)))

# plot - add some empty plots for easier postprocessing in inkscape
{
  pdf(file = "snow_increase_order_MT_compositional.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plot(1,1)
  plot(1,1)
  g <- as.factor(meta_use$treatment)
  r_boxplot <- range(vals, na.rm = T)
  boxplot(vals~g, col=c("gainsboro", "firebrick3", "dodgerblue2"), xlab = "", ylim = r_boxplot, ylab = "Relative abundance (%)", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  dev.off()
}

# snow decrease
snow_decrease_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(3)], sig_values_lmm_best[no_interaction_features,c(3)], fdr_strict_best[no_interaction_features,c(3)]) 
colnames(snow_decrease_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_decrease_res_simple) <- no_interaction_features
de_snow_decrease_simple <- snow_decrease_res_simple[which(snow_decrease_res_simple$fdr<=0.1),]
de_snow_decrease_simple <- de_snow_decrease_simple[order(de_snow_decrease_simple$pvalue),]

# save results into excel
# add worksheet
addWorksheet(wb, "DE_snow_decrease_order_MT")
writeData(wb, sheet = "DE_snow_decrease_order_MT", x = de_snow_decrease_simple, startRow = 1, startCol = 1, rowNames = TRUE)

# plot heatmap
data <- datas_clr$order_mt
data_heat <- data[rownames(de_snow_decrease_simple),]

# define scale and colors
breakslist <- seq(from=-2, to=2, length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))

# add phylum names to plot
tax_mt <- data.frame(tax_table(psq_order_mt_tax), stringsAsFactors = F, check.names = F)
tax_mt <- tax_mt[rownames(data_heat),]

# add GTDB phylum nams
all_phyla <- tax_mt$Phylum
all_phyla <- silva_gtdb_map[all_phyla, 2]
rownames(data_heat) <- paste(rownames(data_heat), all_phyla, sep = "; ")

# order data based on treatment for the heatmap
meta_temp <- meta_use[order(meta_use$treatment),]

# plot only the relevant conditions
meta_temp <- meta_temp[-which(meta_temp$treatment=="+S"),]
data_heat <- data_heat[,rownames(meta_temp)]

# italicize names
ph <- pheatmap(mat = data_heat, cluster_rows = TRUE, cluster_cols = FALSE, scale="row", annotation_col = annotation_col, annotation_colors = annotation_colors, gaps_col = 12, cellwidth = 14, cellheight = 16, breaks=breakslist, color = cols_heat, silent = T, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "")

# extract gtable
gt <- ph$gtable

# identify rowname grobs and set fontface = italic
row_ids <- grep("row_names", gt$layout$name)
for (i in row_ids) {
  gr <- gt$grobs[[i]]
  if (inherits(gr, "text")) {
    gr$gp$fontface <- "italic"
    gt$grobs[[i]] <- gr
  }
}

# Export final PDF
pdf("snow_decrease_order_MT_clr_transformed_only_relevant_samples.pdf", width = 11, height = 4)
grid.newpage()
grid.draw(gt)
dev.off()

# features with significant interactions between the exclusion and snow treatments
interaction_features <- rownames(sig_values_lmm_best)[which(sig_values_lmm_best$best_model=="interaction")]
grazing_treatment_interaction <- data.frame(coefficients_lmm_best[interaction_features,c(4:5)], sig_values_lmm_best[interaction_features,c(4:5)], fdr_strict_best[interaction_features,c(4:5)], stringsAsFactors = F, check.names = F) 
colnames(grazing_treatment_interaction) <- paste(c("coef", "coef", "pvalue", "pvalue", "fdr", "fdr"), colnames(grazing_treatment_interaction))
rownames(grazing_treatment_interaction) <- interaction_features
de_grazing_treatment_interaction <- grazing_treatment_interaction[which(grazing_treatment_interaction$`fdr grazinggrazed:treatment+S`<=0.1 | grazing_treatment_interaction$`fdr grazinggrazed:treatment-S`<=0.1),]

# 3 de features, heatmap, all related to snow decrease and grazing
# save results into excel
# add worksheet
addWorksheet(wb, "DE_grazing_snow_int_order_MT")
writeData(wb, sheet = "DE_grazing_snow_int_order_MT", x = de_grazing_treatment_interaction, startRow = 1, startCol = 1, rowNames = TRUE)

# use boxplots - interactions easier to see - add an extra plot for easier postprocessing in inkscape
tax_mt <- data.frame(tax_table(psq_order_mt_tax), stringsAsFactors = F, check.names = F)
tax_mt <- tax_mt[rownames(de_grazing_treatment_interaction),]

# GTDB annotations
all_phyla <- tax_mt$Phylum
all_phyla <- silva_gtdb_map[all_phyla, 2]

# plot boxplots
{
  pdf(file = "grazing_snow_decrease_interaction_order_MT_compositional_boxplots.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plot(1,1)
  for(ro in 1:nrow(de_grazing_treatment_interaction)){
    vals <- as.numeric(data_order_mt[rownames(de_grazing_treatment_interaction)[ro],])
    
    # add annotations
    main_text <-  paste(all_phyla[ro], rownames(de_grazing_treatment_interaction)[ro], sep = "; ")
    
    # wrap
    main_text <- paste(strwrap(main_text, width = 20), collapse = "\n")
    
    # italize
    main_text <- bquote(italic(.(main_text)))
    
    g <- as.factor(paste(meta_use$treatment, meta_use$grazing, sep = " "))
    g <- factor(g, levels = c("CTL ungrazed", "CTL grazed", "+S ungrazed", "+S grazed", "-S ungrazed", "-S grazed"))
    
    r_boxplot <- range(vals, na.rm = T)
    boxplot(vals~g, col=c("forestgreen", "tan1"), las=2, xlab = "", ylim = r_boxplot, ylab = "Relative abundance (%)", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
  }
  dev.off()
}

# phylum level

# load the lmm results
load(paste(project_root, "/downstream/lmm_analysis/phylum_mt/lmm_analysis_done.RData", sep = ""))

# transform compositional data into percent proportions
data_phylum_mt <- data_phylum_mt * 100

# define the DE lists
# grazing
no_interaction_features <- rownames(sig_values_lmm_best)[-which(sig_values_lmm_best$best_model=="interaction")]
grazing_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(1)], sig_values_lmm_best[no_interaction_features,c(1)], fdr_strict_best[no_interaction_features,c(1)]) 
colnames(grazing_res_simple) <- c("coef", "pvalue", "fdr")
rownames(grazing_res_simple) <- no_interaction_features
de_grazing_simple <- grazing_res_simple[which(grazing_res_simple$fdr<=0.1),]
de_grazing_simple <- de_grazing_simple[order(de_grazing_simple$pvalue),]
# 3 de features

# save results into excel
# add worksheet
addWorksheet(wb, "DE_grazing_phylum_MT")
writeData(wb, sheet = "DE_grazing_phylum_MT", x = de_grazing_simple, startRow = 1, startCol = 1, rowNames = TRUE)

# plot as heatmap
data <- datas_clr$phylum_mt
data_heat <- data[rownames(de_grazing_simple),]

# add gtdb phylum names
all_phyla <- rownames(data_heat)
all_phyla <- silva_gtdb_map[all_phyla, 2]
rownames(data_heat) <- all_phyla

# define scale and colors
breakslist <- seq(from=-2, to=2, length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))

# heatmap
ph <- pheatmap(mat = data_heat, cluster_rows = TRUE, cluster_cols = FALSE, scale="row", annotation_col = annotation_col, annotation_colors = annotation_colors, gaps_col = 18, cellwidth = 14, cellheight = 16, breaks=breakslist, color = cols_heat, silent = T, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "")

# extract gtable
gt <- ph$gtable

# identify rowname grobs and set fontface = italic
row_ids <- grep("row_names", gt$layout$name)
for (i in row_ids) {
  gr <- gt$grobs[[i]]
  if (inherits(gr, "text")) {
    gr$gp$fontface <- "italic"
    gt$grobs[[i]] <- gr
  }
}

# export final PDF
pdf("grazing_phylum_MT_clr_transformed.pdf", width = 11, height = 4)
grid.newpage()
grid.draw(gt)
dev.off()

# snow increase
snow_increase_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(2)], sig_values_lmm_best[no_interaction_features,c(2)], fdr_strict_best[no_interaction_features,c(2)]) 
colnames(snow_increase_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_increase_res_simple) <- no_interaction_features
de_snow_increase_simple <- snow_increase_res_simple[which(snow_increase_res_simple$fdr<=0.1),]
de_snow_increase_simple <- de_snow_increase_simple[order(de_snow_increase_simple$pvalue),]
# no significant findings

# snow decrease
snow_decrease_res_simple <- data.frame(coefficients_lmm_best[no_interaction_features,c(3)], sig_values_lmm_best[no_interaction_features,c(3)], fdr_strict_best[no_interaction_features,c(3)]) 
colnames(snow_decrease_res_simple) <- c("coef", "pvalue", "fdr")
rownames(snow_decrease_res_simple) <- no_interaction_features
de_snow_decrease_simple <- snow_decrease_res_simple[which(snow_decrease_res_simple$fdr<=0.1),]
de_snow_decrease_simple <- de_snow_decrease_simple[order(de_snow_decrease_simple$pvalue),]
# 6 de features

# save results into excel
# add worksheet
addWorksheet(wb, "DE_snow_decrease_phylum_MT")
writeData(wb, sheet = "DE_snow_decrease_phylum_MT", x = de_snow_decrease_simple, startRow = 1, startCol = 1, rowNames = TRUE)

# plot as heatmap
data <- datas_clr$phylum_mt
data_heat <- data[rownames(de_snow_decrease_simple),]

# order data based on treatment for the heatmap
meta_temp <- meta_use[order(meta_use$treatment),]

# plot only the relevant conditions
meta_temp <- meta_temp[-which(meta_temp$treatment=="+S"),]
data_heat <- data_heat[,rownames(meta_temp)]

# add GTDB annotations
all_phyla <- rownames(data_heat)
all_phyla <- silva_gtdb_map[all_phyla, 2]
rownames(data_heat) <- all_phyla

ph <- pheatmap(mat = data_heat, cluster_rows = TRUE, cluster_cols = FALSE, scale="row", annotation_col = annotation_col, annotation_colors = annotation_colors, gaps_col = 12, cellwidth = 14, cellheight = 16, breaks=breakslist, color = cols_heat, silent = T, fontsize_row = 8, fontsize_col = 8, filename = NA, main = "")

# extract gtable
gt <- ph$gtable

# identify rowname grobs and set fontface = italic
row_ids <- grep("row_names", gt$layout$name)
for (i in row_ids) {
  gr <- gt$grobs[[i]]
  if (inherits(gr, "text")) {
    gr$gp$fontface <- "italic"
    gt$grobs[[i]] <- gr
  }
}

# export final PDF
pdf("snow_decrease_phylum_MT_clr_transformed_only_relevant_samples.pdf", width = 8, height = 3)
grid.newpage()
grid.draw(gt)
dev.off()

# features with significant interactions between the exclusion and snow treatments
interaction_features <- rownames(sig_values_lmm_best)[which(sig_values_lmm_best$best_model=="interaction")]
grazing_treatment_interaction <- data.frame(coefficients_lmm_best[interaction_features,c(4:5)], sig_values_lmm_best[interaction_features,c(4:5)], fdr_strict_best[interaction_features,c(4:5)], stringsAsFactors = F, check.names = F) 
colnames(grazing_treatment_interaction) <- paste(c("coef", "coef", "pvalue", "pvalue", "fdr", "fdr"), colnames(grazing_treatment_interaction))
rownames(grazing_treatment_interaction) <- interaction_features
de_grazing_treatment_interaction <- grazing_treatment_interaction[which(grazing_treatment_interaction$`fdr grazinggrazed:treatment+S`<=0.1 | grazing_treatment_interaction$`fdr grazinggrazed:treatment-S`<=0.1),]
# one de feature

# save results into excel
# add worksheet
addWorksheet(wb, "DE_grazing_snow_int_phylum_MT")
writeData(wb, sheet = "DE_grazing_snow_int_phylum_MT", x = de_grazing_treatment_interaction, startRow = 1, startCol = 1, rowNames = TRUE)

# one feature, use boxplot - interactions more easy to see - add extra empty plots for easier postprocessing in inkscape

# map to silva phylum
all_phyla <- rownames(de_grazing_treatment_interaction)
all_phyla <- silva_gtdb_map[all_phyla, 2]

# plot boxplot
{
  pdf(file = "grazing_snow_decrease_interaction_phylum_MT_compositional_boxplots.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plot(1,1)
  for(ro in 1:nrow(de_grazing_treatment_interaction)){
    vals <- as.numeric(data_phylum_mt[rownames(de_grazing_treatment_interaction)[ro],])
    
    # add annotations
    main_text <-  rownames(de_grazing_treatment_interaction)[ro]
    
    # italize
    main_text <- bquote(italic(.(main_text)))
    
    g <- as.factor(paste(meta_use$treatment, meta_use$grazing, sep = " "))
    g <- factor(g, levels = c("CTL ungrazed", "CTL grazed", "+S ungrazed", "+S grazed", "-S ungrazed", "-S grazed"))
    
    r_boxplot <- range(vals, na.rm = T)
    boxplot(vals~g, col=c("forestgreen", "tan1"), las=2, xlab = "", ylim = r_boxplot, ylab = "Relative abundance (%)", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
  }
  dev.off()
}

# save compiled excel workbook
saveWorkbook(wb, "figure3_DA_DE_taxa.xlsx", overwrite = TRUE)

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures