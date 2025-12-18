# a small script to plot the abundance of MAGs in different samples

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(phyloseq)
library(microbiome)
library(grid)
library(RColorBrewer)
library(pheatmap)

# load taxonomy - use the modified taxonomy generated in the script plotting the most abundant taxonomic annotations for the MAGs
load(paste(project_root, "/downstream/manuscript_figures/figure6/modified_gtdb_tax_silva.RData", sep = ""))

# load MAG abundance data - coverM
setwd(paste(project_root, "/metagenomics/mag_based/final_mags/coverm", sep = ""))

# load coverM results
coverm_mg <- read.csv("coverM_MG_final_metrics.tsv", sep = "\t", stringsAsFactors = F, check.names = F)
coverm_mg <- coverm_mg[-1,]
rownames(coverm_mg) <- coverm_mg$Genome

# use relative abundance
coverm_mg <- coverm_mg[,grep("Relative Abundance", colnames(coverm_mg))]

# parse the sample names
colnames(coverm_mg) <- gsub("_L004_R1_trimmed.fastq Relative Abundance", "", colnames(coverm_mg))
colnames(coverm_mg) <- gsub("[ (%)]", "", colnames(coverm_mg))
colnames(coverm_mg) <- paste("P",unlist(lapply(colnames(coverm_mg), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

# load metadata and organize datas acccordingly
load(paste(project_root, "/metadata/Study_Metadata.RData", sep = ""))
rownames(metadata) <- metadata$`Short code`
coverm_mg <- coverm_mg[,rownames(metadata)]

# add MAG taxonomy as rownames: MAG_number;Phyla(GTDB Phyla);Order;Genus
gtdbtk_report <- gtdbtk_report[rownames(coverm_mg),]
m_num <- seq(1:nrow(gtdbtk_report))
r_names <- paste("MAG",m_num, "; ",gtdbtk_report$Phylum, "; ", gtdbtk_report$Order, "; ", gtdbtk_report$Genus, sep = "")

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# annotations for plotting 
annotation_col <- meta_use
for(j in 1:ncol(annotation_col)){annotation_col[,j] <- as.character(annotation_col[,j])}

# use the same colors as for the other heatmaps in the manuscript
annotation_colors <- list(
  grazing = c("grazed"="tan1", "ungrazed"="forestgreen"),
  veg_clusters = c("t_ces"="darkolivegreen2", "c_ros"="purple", "c_cho"="deeppink"),
  treatment = c("-S"="dodgerblue2", "CTL"="gainsboro", "+S"="firebrick3")
)

# plot
# define range for colors
r <- as.numeric(unlist(coverm_mg))
r <- c(quantile(r, 0.05, na.rm = T), quantile(r, 0.95, na.rm = T))

# define colors
breakslist <- seq(from=r[1], to=r[2], length.out=11)
cols_heat <- rev(brewer.pal(11,"RdBu"))

data_heat <- coverm_mg
rownames(data_heat) <- r_names

# plot heatmap 
ph <- pheatmap(mat = data_heat, na_col = "darkgrey", clustering_distance_rows = "correlation", fontsize = 12, cluster_rows = T, breaks = breakslist, color = cols_heat, cluster_cols = F, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="none", cellwidth = 14, cellheight = 14, silent = T, fontsize_row = 14, fontsize_col = 14, filename = NA)

# italicize MAG names
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

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure6", sep = ""))

# export final pdf
pdf("MAG_MG_Abundance.pdf", width = 16, height = 24)
grid.newpage()
grid.draw(gt)
dev.off()

# plot number of marker gene copies for the same MAGs
# load this info
load(paste(project_root, "/metagenomics/mag_based/final_mags/Final_MAG_Taxonomy_Mark_Gene.RData", sep = ""))

# get the genes of interst
int_genes <- c("NirS", "NirK", "NorB", "NosZ", "PmoA", "McrA")
int_gene_copes <- t(markergene_copies_clusters[int_genes, ])

# plot

# define range for colors
r <- as.numeric(unlist(int_gene_copes))
r <- range(r)

# define colors
cols_heat <- brewer.pal(r[2]+1,"Reds")

# copy numbers going from 0 to 4
breakslist <- seq(from=0, to=4, length.out=6)

# plot heatmap
int_gene_copes <- int_gene_copes[rownames(coverm_mg),]
data_heat <- int_gene_copes
rownames(data_heat) <- r_names

# organize according to the clustering in the abundance heatmap
c_names <- colnames(data_heat)
ann <- data.frame(Gene = c_names)
rownames(ann) <- c_names
data_heat <- data_heat[ph$tree_row$order,]

# plot heatmap into object
ph <- pheatmap(mat = data_heat, na_col = "darkgrey", fontsize = 12, gaps_col = 4, cluster_rows = F, breaks = breakslist, color = cols_heat, cluster_cols = F, scale="none", annotation_col = ann, cellwidth = 14, cellheight = 14, silent = T, fontsize_row = 14, show_colnames = T, fontsize_col = 14, filename = NA)

# italicize MAG names
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

# export final pdf
pdf("N_CH4_marker_genes_nr_copies_MAGs_MG.pdf", width = 10, height = 24)
grid.newpage()
grid.draw(gt)
dev.off()

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures
