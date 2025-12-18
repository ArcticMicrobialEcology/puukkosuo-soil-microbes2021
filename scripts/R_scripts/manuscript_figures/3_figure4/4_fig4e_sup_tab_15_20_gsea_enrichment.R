# a small script to plot the results of the gene set enrichment analysis (GSEA) for the KEGG modules

# load libraries
library(openxlsx)

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load and plot fgsea results
# metagenomics
load(paste(project_root, "/downstream/lmm_analysis/kegg_gsea_enrichment/kegg_mg/fgsea_results.RData", sep = ""))

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure4", sep = ""))

# grazing / exclusion treatment
temp_fgsea <- fgsea_grazed
temp_fgsea <- temp_fgsea[order(as.numeric(temp_fgsea$pval), decreasing = F),]

# pick only the significant results, fdr <= 0.1
sig_fgsea <- temp_fgsea[which(temp_fgsea$padj<=0.1),]

# save table into excel
# create a workbook
wb <- createWorkbook()

# save table into excel
# add worksheet
addWorksheet(wb, "GSEA_grazing_MG")
writeData(wb, sheet = "GSEA_grazing_MG", x = sig_fgsea, startRow = 1, startCol = 1, rowNames = TRUE)

# plot as barplot
# define colors
cols_use <- rep("firebrick1", nrow(sig_fgsea))
cols_use[which(sig_fgsea$NES<0)] <- "blue"

{
  pdf(file ="grazing_KEGG_MG.pdf", width = 10, height = 3, onefile = T)
  par(mar=c(6.1, 35.1, 4.1, 2.1))
  barplot(sig_fgsea$padj, horiz = TRUE, names.arg = gsub("_", " ", sig_fgsea$pathway), las=1, xlim = c(0, 0.12), col = cols_use, xlab = "Adjusted P-value")
  abline(v = 0.1, lwd=2)
  abline(v = 0.05, lwd=2, lty=2)
  dev.off()
}

# snow addition
temp_fgsea <- fgsea_snow_addition
temp_fgsea <- temp_fgsea[order(as.numeric(temp_fgsea$pval), decreasing = F),]

# pick only the significant results, fdr <= 0.1
sig_fgsea <- temp_fgsea[which(temp_fgsea$padj<=0.1),]

# save table into excel
# add worksheet
addWorksheet(wb, "GSEA_snow_addition_MG")
writeData(wb, sheet = "GSEA_snow_addition_MG", x = sig_fgsea, startRow = 1, startCol = 1, rowNames = TRUE)

# plot as barplot
# define colors
cols_use <- rep("firebrick1", nrow(sig_fgsea))
cols_use[which(sig_fgsea$NES<0)] <- "blue"

{
  pdf(file ="snow_addition_KEGG_MG.pdf", width = 10, height = 5, onefile = T)
  par(mar=c(6.1, 35.1, 4.1, 2.1))
  barplot(sig_fgsea$padj, horiz = TRUE, names.arg = gsub("_", " ", sig_fgsea$pathway), las=1, xlim = c(0, 0.12), col = cols_use, xlab = "Adjusted P-value")
  abline(v = 0.1, lwd=2)
  abline(v = 0.05, lwd=2, lty=2)
  dev.off()
}

# snow removal
temp_fgsea <- fgsea_snow_decrease
temp_fgsea <- temp_fgsea[order(as.numeric(temp_fgsea$pval), decreasing = F),]

# pick only the significant results, fdr <= 0.1
sig_fgsea <- temp_fgsea[which(temp_fgsea$padj<=0.1),]

# save table into excel
# add worksheet
addWorksheet(wb, "GSEA_snow_removal_MG")
writeData(wb, sheet = "GSEA_snow_removal_MG", x = sig_fgsea, startRow = 1, startCol = 1, rowNames = TRUE)

# plot as barplot
# define colors
cols_use <- rep("firebrick1", nrow(sig_fgsea))
cols_use[which(sig_fgsea$NES<0)] <- "blue"

{
  pdf(file ="snow_removal_KEGG_MG.pdf", width = 10, height = 3.5, onefile = T)
  par(mar=c(6.1, 35.1, 4.1, 2.1))
  barplot(sig_fgsea$padj, horiz = TRUE, names.arg = gsub("_", " ", sig_fgsea$pathway), las=1, xlim = c(0, 0.12), col = cols_use, xlab = "Adjusted P-value")
  abline(v = 0.1, lwd=2)
  abline(v = 0.05, lwd=2, lty=2)
  dev.off()
}

# metatranscriptomics
load(paste(project_root, "/downstream/lmm_analysis/kegg_gsea_enrichment/kegg_mt/fgsea_results.RData", sep = ""))

# grazing / exclusion treatment
temp_fgsea <- fgsea_grazed
temp_fgsea <- temp_fgsea[order(as.numeric(temp_fgsea$pval), decreasing = F),]

# pick only the significant results, fdr <= 0.1
sig_fgsea <- temp_fgsea[which(temp_fgsea$padj<=0.1),]

# save table into excel
# add worksheet
addWorksheet(wb, "GSEA_grazing_MT")
writeData(wb, sheet = "GSEA_grazing_MT", x = sig_fgsea, startRow = 1, startCol = 1, rowNames = TRUE)

# plot as barplot
# define colors
cols_use <- rep("firebrick1", nrow(sig_fgsea))
cols_use[which(sig_fgsea$NES<0)] <- "blue"

{
  pdf(file ="grazing_KEGG_MT.pdf", width = 10, height = 5, onefile = T)
  par(mar=c(6.1, 35.1, 4.1, 2.1))
  barplot(sig_fgsea$padj, horiz = TRUE, names.arg = gsub("_", " ", sig_fgsea$pathway), las=1, xlim = c(0, 0.12), col = cols_use, xlab = "Adjusted P-value")
  abline(v = 0.1, lwd=2)
  abline(v = 0.05, lwd=2, lty=2)
  dev.off()
}

# snow addition
temp_fgsea <- fgsea_snow_addition
temp_fgsea <- temp_fgsea[order(as.numeric(temp_fgsea$pval), decreasing = F),]

# pick only the significant results, fdr <= 0.1
sig_fgsea <- temp_fgsea[which(temp_fgsea$padj<=0.1),]

# save table into excel
# add worksheet
addWorksheet(wb, "GSEA_snow_addition_MT")
writeData(wb, sheet = "GSEA_snow_addition_MT", x = sig_fgsea, startRow = 1, startCol = 1, rowNames = TRUE)

# plot as barplot
# define colors
cols_use <- rep("firebrick1", nrow(sig_fgsea))
cols_use[which(sig_fgsea$NES<0)] <- "blue"

{
  pdf(file ="snow_addition_KEGG_MT.pdf", width = 10, height = 5, onefile = T)
  par(mar=c(6.1, 35.1, 4.1, 2.1))
  barplot(sig_fgsea$padj, horiz = TRUE, names.arg = gsub("_", " ", sig_fgsea$pathway), las=1, xlim = c(0, 0.12), col = cols_use, xlab = "Adjusted P-value")
  abline(v = 0.1, lwd=2)
  abline(v = 0.05, lwd=2, lty=2)
  dev.off()
}

# snow removal
temp_fgsea <- fgsea_snow_decrease
temp_fgsea <- temp_fgsea[order(as.numeric(temp_fgsea$pval), decreasing = F),]

# pick only the significant results, fdr <= 0.1
sig_fgsea <- temp_fgsea[which(temp_fgsea$padj<=0.1),]

# save table into excel
# add worksheet
addWorksheet(wb, "GSEA_snow_removal_MT")
writeData(wb, sheet = "GSEA_snow_removal_MT", x = sig_fgsea, startRow = 1, startCol = 1, rowNames = TRUE)

# plot as barplot
# define colors
cols_use <- rep("firebrick1", nrow(sig_fgsea))
cols_use[which(sig_fgsea$NES<0)] <- "blue"

{
  pdf(file ="snow_removal_KEGG_MT.pdf", width = 10, height = 2.5, onefile = T)
  par(mar=c(6.1, 35.1, 4.1, 2.1))
  barplot(sig_fgsea$padj, horiz = TRUE, names.arg = gsub("_", " ", sig_fgsea$pathway), las=1, xlim = c(0, 0.12), col = cols_use, xlab = "Adjusted P-value")
  abline(v = 0.1, lwd=2)
  abline(v = 0.05, lwd=2, lty=2)
  dev.off()
}

# save compiled excel workbook
saveWorkbook(wb, "figure4_KEGG_GSEA_SIG.xlsx", overwrite = TRUE)

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures