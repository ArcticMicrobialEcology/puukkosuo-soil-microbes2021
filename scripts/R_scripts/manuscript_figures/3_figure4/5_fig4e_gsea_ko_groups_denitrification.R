# a small script to plot the results of the gene set enrichment analysis (GSEA) 
# in the module M00529_Denitrification, nitrate => nitrogen on a ko group level

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load data
# metagenomics
load(paste(project_root, "/metagenomics/kegg_diamond/Matrices_For_Downstream.RData", sep = ""))

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

# load the metagenomics fgsea results
load(paste(project_root, "/downstream/lmm_analysis/kegg_gsea_enrichment/kegg_mg/fgsea_results.RData", sep = ""))

# grazing / exclusion treatment fgsea results
temp_fgsea <- fgsea_grazed
temp_fgsea <- temp_fgsea[order(as.numeric(temp_fgsea$pval), decreasing = F),]

# module M00529_Denitrification, nitrate => nitrogen is the most enriched module, plot only that one
sig_fgsea <- temp_fgsea[1,]

# use tpm data
data <- ko_tpm_data

# set colors for plotting
col_grazed <- "tan1"
col_ungrazed <- "forestgreen"

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure4", sep = ""))

{
  pdf(file = "MG_gsea_M00529_ko_groups.pdf", width = 11.7, height = 8.3 , onefile = T)
  par(mar=c(15.1, 5.1, 4.1, 2.1))
  
  # get the ko groups
  ko_int <- gsub(" ", "", unlist(strsplit(x = sig_fgsea$leadingEdge[1], ";")))
  
  # get the related data
  data_temp <- data[ko_int,]
  
  # scale the abundance of each ko group to put in the same plot
  # scale each ko’s tpm across samples (mean 0, sd 1) for plotting
  data_temp <- t(scale(t(data_temp)))
  annot_plot <- ko_annotations[rownames(data_temp),]
  annot_plot <- unlist(lapply(strsplit(x = annot_plot, split = ";"), function(x) x[1]))
  rownames(data_temp) <- paste(rownames(data_temp), annot_plot, sep = "; ")
  
  # create an empty plot
  plot(1,
       xlim = c(0, (nrow(data_temp)*3)+2),
       xaxt = "n",
       ylim = c(-4,4),
       xlab = "", 
       ylab = "Scaled copies per million",
       main = gsub("_", " ", sig_fgsea$pathway[1]),
       bty = "n",
       cex.axis = 1.5,
       cex.lab = 1.5,
       type = "n") 
  
  ind1 <- 1
  ind2 <- 2
  ind_lab <- 1.5
  x_axt_tics <- numeric(0)
  
  # add all the ko groups
  for(feat in 1:nrow(data_temp)){
    vals <- data_temp[feat,]
    g <- meta_use$grazing
    boxplot(vals~g, col=c(col_ungrazed, col_grazed), add = T, at=c(ind1, ind2), xlab = "", xaxt="n", ylab = "", yaxt="n", frame=FALSE)
    
    ind1 <- ind1 + 3
    ind2 <- ind2 + 3
    x_axt_tics <- c(x_axt_tics, ind_lab)
    ind_lab <- ind_lab +3
  }
  axis(side = 1, at = x_axt_tics, labels = rownames(data_temp), las=2)
  dev.off()
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures


