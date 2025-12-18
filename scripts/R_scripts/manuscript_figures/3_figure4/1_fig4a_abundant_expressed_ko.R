# a small script to plot the most abundant and expressed KO groups in the MG and MT datas

# load libraries
library(grid)
library(gridExtra)

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load datas 
# metagenomics
load(paste(project_root, "/metagenomics/kegg_diamond/Matrices_For_Downstream.RData", sep = ""))
data_mg <- ko_tpm_data

# metatranscriptomics
load(paste(project_root, "/metatranscriptomics/kegg_diamond/Matrices_For_Downstream.RData", sep = ""))
data_mt <- ko_tpm_data

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

# define some divergent colors for plotting
annotation_colors <- c("purple", "cornflowerblue","mistyrose2","darkgrey","deeppink","brown","violet","firebrick1","chocolate3","forestgreen","darkslateblue","darkolivegreen",
                       "chocolate1", "cadetblue1", "darkgoldenrod3", "blue", "darkolivegreen3", "yellow", "burlywood2", "black", "green")

# get 10 most expressed features over all samples on average
# transform into closed proportions
c_sums <- colSums(data_mg)
data_prop_mg <- data_mg
for(j in 1:ncol(data_prop_mg)){data_prop_mg[,j] <- data_prop_mg[,j] / c_sums[j]}
data_prop_mg <- data_prop_mg * 100

c_sums <- colSums(data_mt)
data_prop_mt <- data_mt
for(j in 1:ncol(data_prop_mt)){data_prop_mt[,j] <- data_prop_mt[,j] / c_sums[j]}
data_prop_mt <- data_prop_mt * 100

# feature means in proportional data
feature_means_prop_mg <- rowMeans(data_prop_mg, na.rm = T)
feature_means_prop_mg <- feature_means_prop_mg[order(feature_means_prop_mg, decreasing = T)]

feature_means_prop_mt <- rowMeans(data_prop_mt, na.rm = T)
feature_means_prop_mt <- feature_means_prop_mt[order(feature_means_prop_mt, decreasing = T)]

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure4", sep = ""))

# use the same range for both MG and MT
r <- c(0, 14.5)

# metagenomics
{
  pdf(file = "Most_expressed_features_samplewise_on average_MG_prop_same_range.pdf", width = 8.3, height = 11.7, onefile = T)
  
  # desired bar width in inches
  desired_bar_width_inches <- 1
  
  # total page width in inches
  total_page_width_inches <- 8.3
  
  # calculate the left and right margins to center the bar
  margin_inches <- (total_page_width_inches - desired_bar_width_inches) / 2
  
  # set the plot dimensions in inches
  options(repr.plot.width = total_page_width_inches, repr.plot.height = 4, repr.plot.units = "in") 
  
  # calculate relative bar width and margins
  bar_width_relative <- desired_bar_width_inches / total_page_width_inches
  margin_relative <- margin_inches / total_page_width_inches
  
  # plot
  barplot(as.matrix(feature_means_prop_mg[1:10]), main = "" , width = bar_width_relative, col=annotation_colors[1:10],space = c(margin_relative, 0), xlim = c(0,1), ylim = r, border="white", las=2, xlab="Sample", ylab = "Proportion (%)")
  
  # reset graphics options
  options(repr.plot.width = NULL, repr.plot.height = NULL, repr.plot.units = NULL)
  
  # add legend to the next page
  grid.newpage()
  legend("top", legend = names(feature_means_prop_mg)[1:10], fill = annotation_colors[1:10], ncol = 4, bty="n")
  
  annots_plot <- character(10)
  for(ro in 1:10){annots_plot[ro] <-  paste(c(names(feature_means_prop_mg)[ro], as.character(annot_use[names(feature_means_prop_mg)[ro],])), collapse = "; ")}
  
  df<-data.frame(text = annots_plot)
  d  <- sapply(lapply(df$text, strwrap, width=250), paste, collapse="\n")
  grid.draw(tableGrob(d,theme=ttheme_minimal(base_size = 6) ))
  dev.off()
}

# metatranscriptomics
{
  pdf(file = "Most_expressed_features_samplewise_on average_MT_prop_same_range.pdf", width = 8.3, height = 11.7, onefile = T)
  
  # desired bar width in inches
  desired_bar_width_inches <- 1
  
  # total page width in inches
  total_page_width_inches <- 8.3
  
  # calculate the left and right margins to center the bar
  margin_inches <- (total_page_width_inches - desired_bar_width_inches) / 2
  
  # set the plot dimensions in inches
  options(repr.plot.width = total_page_width_inches, repr.plot.height = 4, repr.plot.units = "in") 
  
  # calculate relative bar width and margins
  bar_width_relative <- desired_bar_width_inches / total_page_width_inches
  margin_relative <- margin_inches / total_page_width_inches
  
  # plot
  barplot(as.matrix(feature_means_prop_mt[1:10]), main = "" , width = bar_width_relative, ylim = r, col=annotation_colors[11:20],space = c(margin_relative, 0), xlim = c(0,1), border="white", las=2, xlab="Sample", ylab = "Proportion (%)")
  
  # reset graphics options
  options(repr.plot.width = NULL, repr.plot.height = NULL, repr.plot.units = NULL)
  
  # add legend to the next page
  grid.newpage()
  legend("top", legend = names(feature_means_prop_mt)[1:10], fill = annotation_colors[11:20], ncol = 4, bty="n")
  
  annots_plot <- character(10)
  for(ro in 1:10){annots_plot[ro] <-  paste(c(names(feature_means_prop_mt)[ro], as.character(annot_use[names(feature_means_prop_mt)[ro],])), collapse = "; ")}
  
  df<-data.frame(text = annots_plot)
  d  <- sapply(lapply(df$text, strwrap, width=250), paste, collapse="\n")
  grid.draw(tableGrob(d,theme=ttheme_minimal(base_size = 6) ))
  dev.off()
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures
