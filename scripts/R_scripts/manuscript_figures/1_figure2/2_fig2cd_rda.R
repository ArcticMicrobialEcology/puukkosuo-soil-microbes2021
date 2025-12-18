# a small script to plot RDA ordinations for various datas

# define some needed directories
args <- character(4)
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# kegg KO MG data - processed filtered by the processing scripts
args[1] <- paste(project_root, "/metagenomics/kegg_diamond/Matrices_For_Downstream.RData", sep = "")

# phyloflash MG taxonomy data - processed filtered by the processing scripts
args[2] <- paste(project_root, "/metagenomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = "")

# kegg KO MT data - processed filtered by the processing scripts
args[3] <- paste(project_root, "/metatranscriptomics/kegg_diamond/Matrices_For_Downstream.RData", sep = "")

# phyloflash MT taxonomy data - processed filtered by the processing scripts
args[4] <- paste(project_root, "/metatranscriptomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = "")

# define plotting directory - where the plots will be produced - needs to exist

plotting_dir <- paste(project_root, "downstream/manuscript_figures/figure2", sep = "")

# load libraries
library(vegan)
library(phyloseq)

# define some custom functions to be used
# a function to adjust color transparency
makeTransparent = function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}

# a function for RDA with arrows from regression coefficients and varying point shape and size for vegetation and snow treatment
plotRDAreg <- function(data_for_pca, metadata, colors, legend_graz_loc="topleft", legend_veg_loc="topright", legend_treat_loc="bottomleft"){ # color for grazed first, then color for ungrazed
  
  # pick the relevant metadata
  meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
  colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")
  
  # set levels for metadata
  meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
  meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
  meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))
  
  # perform rda
  set.seed(1)
  rda_data <- rda(t(data_for_pca) ~ ., meta_use)
  
  # define colors
  cols <- rep(colors[1], nrow(meta_use))
  cols[which(meta_use$grazing=="ungrazed")] <- colors[2]
  
  # define point shapes
  pchs <- rep(19, nrow(meta_use))
  pchs[which(meta_use$veg_clusters=="c_cho")] <- 15
  pchs[which(meta_use$veg_clusters=="c_ros")] <- 17
  
  # define point sizes - snow treatment
  cexes <- rep(2, nrow(meta_use))
  cexes[which(meta_use$treatment=="-S")] <- 1
  cexes[which(meta_use$treatment=="+S")] <- 3
  
  # define ranges for the x and y x axes, not always very nice by default
  xlim_rda <- range(scores(rda_data)$sites[,1])
  ylim_rda <- range(scores(rda_data)$sites[,2])
  
  par(mar=c(5.1, 5.1, 4.1, 2.1))
  plot(rda_data, type = "none", main="Redundancy analysis", cex.axis=1.5, cex.lab=1.5, bty="n", xlim=xlim_rda, ylim=ylim_rda) # empty plot
  
  # add points
  points(rda_data, display = "sites", pch = pchs, col = cols, cex=cexes)
  
  # add regression arrows
  text(rda_data, display = c("bp"), col="black", labels =  c("grazed", "S+", "S-", "C_cho", "C_ros"), cex=1.5)
  
  # add legends
  legend(legend_graz_loc, legend = c("grazed", "ungrazed"), fill = colors, cex = 1.5)
  legend(legend_veg_loc, legend = c("t_ces", "c_cho", "c_ros") , col="black", pch = c(19,15,17), cex=1.5)
  legend(legend_treat_loc, legend = c("-S", "CTL", "+S"), col="black", pch = 19, pt.cex=c(1,2,3))
  
}

# define custom color
trans_grey <- makeTransparent("grey", alpha = 0.8)

## metagenomics
# kegg ko data 
load(args[1])

# log2 transform the data for RDA using Euclidean distance (in effect PCA)
data_for_pca <- log2(ko_tpm_data+1)

# change to plotting dir
setwd(plotting_dir)

# perform and plot RDA
{
  pdf(file = "kegg_MG_RDA.pdf", width = 10, height = 10)
  plotRDAreg(data_for_pca = data_for_pca, metadata = metadata, colors = c("tan1", "forestgreen"), legend_graz_loc = "topright", legend_veg_loc = "bottomright")
  dev.off()
}

# taxonomy data
load(args[2])

# use the clr transformed non-aggregated data for RDA
data_for_pca <- otu_table_clr_filt

# perform and plot RDA
{
  pdf(file = "tax_MG_RDA.pdf", width = 10, height = 10)
  plotRDAreg(data_for_pca = data_for_pca, metadata = metadata, colors = c("tan1", "forestgreen"), legend_graz_loc = "topright", legend_veg_loc = "bottomright")
  dev.off()
}

## metatranscriptomics
# kegg 
load(args[3])

# log2 transform the data for RDA using Euclidean distance
data_for_pca <- log2(ko_tpm_data+1)

# perform and plot RDA
{
  pdf(file = "kegg_MT_RDA.pdf", width = 10, height = 10)
  plotRDAreg(data_for_pca = data_for_pca, metadata = metadata, colors = c("tan1", "forestgreen"), legend_graz_loc = "topright", legend_veg_loc = "bottomright", legend_treat_loc = "topleft")
  dev.off()
}

# taxonomy data
load(args[4])

data_for_pca <- otu_table_clr_filt

# perform and plot RDA
{
  pdf(file = "tax_MT_RDA.pdf", width = 10, height = 10)
  plotRDAreg(data_for_pca = data_for_pca, metadata = metadata, colors = c("tan1", "forestgreen"), legend_graz_loc = "bottomleft", legend_veg_loc = "bottomright")
  dev.off()
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures