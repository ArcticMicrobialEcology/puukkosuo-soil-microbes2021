# a small script to plot the expression of selected nitrogen and methane related metabolic marker genes over the discovered MAGs

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(phyloseq)
library(microbiome)
library(grid)
library(RColorBrewer)
library(lmerTest)
library(pheatmap)

# load the summarized metatranscriptomic count datas for MAG genes
load(paste(project_root, "/metagenomics/mag_based/final_mags/MT_Summarized_Count_Datas_For_All_Bins.RData", sep = ""))

# get all the unique marker genes present in the MAG datas
all_mag_marker_genes <- character()
for(i in 1:length(bin_mt_datas)){
  temp <- bin_mt_datas[[i]]
  temp <- temp$tpm$metmarkdb
  all_mag_marker_genes <- c(all_mag_marker_genes, rownames(temp))
}
all_mag_marker_genes <- unique(all_mag_marker_genes)

# prepare a data frame for those marker genes
all_mag_marker_gene_data <- data.frame(matrix(nrow = length(all_mag_marker_genes), ncol = ncol(temp)))
rownames(all_mag_marker_gene_data) <- all_mag_marker_genes
colnames(all_mag_marker_gene_data) <- colnames(temp)
all_mag_marker_gene_data[is.na(all_mag_marker_gene_data)] <- 0

# sum the tpm values of the MAGs for the marker genes together to look at marker genes over all MAGs
for(i in 1:length(bin_mt_datas)){
  temp <- bin_mt_datas[[i]]
  temp <- temp$tpm$metmarkdb
  temp <- temp[,colnames(all_mag_marker_gene_data)]
  all_mag_marker_gene_data[rownames(temp),] <- all_mag_marker_gene_data[rownames(temp),]+temp
}

# parse sample names
colnames(all_mag_marker_gene_data) <- paste("P",unlist(lapply(colnames(all_mag_marker_gene_data), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")

# load metadata and organize datas acccordingly
load(paste(project_root, "/metadata/Study_Metadata.RData", sep = ""))
rownames(metadata) <- metadata$`Short code`
all_mag_marker_gene_data <- all_mag_marker_gene_data[,rownames(metadata)]

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure6", sep = ""))

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# plot the nitrogen and methane related marker genes
# define nitrogen and methane related genes
nitrogen_genes <- c("NirS", "NirK", "NorB", "NosZ","NifH", "NapA", "NrfA", "NarG")
methane_genes <- c("PmoA", "McrA")

# start plotting
data <- all_mag_marker_gene_data

# nitrogen genes
nit_genes <- nitrogen_genes[which(nitrogen_genes%in%rownames(data))]
nit_data <- data[nit_genes,]

# get significance values for the exclusion / grazing treatment
lmm_nitrogen <- list()
sig_val_nit <- data.frame(matrix(nrow = nrow(nit_data), ncol = 2))

for(i in 1:nrow(nit_data)){
  temp_data <- cbind(log2(as.numeric(nit_data[i,])+1), meta_use)
  colnames(temp_data)[1] <- "expression"
  
  # use the simple model which has been the best model for other analysis
  set.seed(1)
  mod_simple <- lmer(expression ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
  
  sum_simple <- summary(mod_simple)
  sig_val_nit[i,1] <- sum_simple$coefficients[2,1] # coefficient
  sig_val_nit[i,2] <- sum_simple$coefficients[2,5] # p-value
  
  # save the models
  lmm_nitrogen[[i]] <- mod_simple
}
rownames(sig_val_nit) <- rownames(nit_data)
colnames(sig_val_nit) <- c("coefficient", "pval")
names(lmm_nitrogen) <- rownames(nit_data)

# gene names for plotting
plot_names <- c("nirS","nirK","norB", "nosZ", "nifH", "napA", "nrfA","narG")

# plot several in one page
{
  pdf(file = "grazing_tpm_values_boxplot_all_one_page_nitrogen.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  for(ro in 1:nrow(nit_data)){
    vals <- as.numeric(nit_data[ro,])
    r_boxplot <- range(vals, na.rm = T)
    main_text <- plot_names[ro]
    
    # italize
    main_text <- bquote(italic(.(main_text)))
    
    g <- as.factor(meta_use$grazing)
    boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, outline = F, ylab = "Copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, frame = FALSE)
    
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
    
    # add sig values
    sig_values_text <- paste("P-value:",round(sig_val_nit[rownames(nit_data)[ro],2],3))
    mtext(sig_values_text, cex=0.8)
  }
  dev.off()
}

# methane
meth_genes <- methane_genes[which(methane_genes%in%rownames(data))]
meth_data <- data[meth_genes,]

# get significance values for grazing
lmm_methane <- list()
sig_val_meth <- data.frame(matrix(nrow = nrow(meth_data), ncol = 2))

for(i in 1:nrow(meth_data)){
  temp_data <- cbind(log2(as.numeric(meth_data[i,])+1), meta_use)
  colnames(temp_data)[1] <- "expression"
  
  # use the simple model which has been the best model for other analysis
  set.seed(1)
  mod_simple <- lmer(expression ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
  
  # only use the simple model here
  sum_simple <- summary(mod_simple)
  sig_val_meth[i,1] <- sum_simple$coefficients[2,1] # coefficient
  sig_val_meth[i,2] <- sum_simple$coefficients[2,5] # p-value
  
  # save the models
  lmm_methane[[i]] <- mod_simple
  
}
rownames(sig_val_meth) <- rownames(meth_data)
colnames(sig_val_meth) <- c("coefficient", "pval")
names(lmm_methane) <- rownames(meth_data)

# gene names properly for plots - mcrA not present
plot_names <- c("pmoA")

# plot similarly
{
  pdf(file = "grazing_tpm_values_boxplot_all_one_page_methane.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  for(ro in 1:nrow(meth_data)){
    vals <- as.numeric(meth_data[ro,])
    r_boxplot <- range(vals, na.rm = T)
    main_text <- plot_names[ro]
    
    # italize
    main_text <- bquote(italic(.(main_text)))
    
    g <- as.factor(meta_use$grazing)
    boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, outline = F, ylab = "Copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, frame = FALSE)
    
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
    
    # add sig values
    sig_values_text <- paste("P-value:",round(sig_val_meth[rownames(meth_data)[ro],2],3))
    mtext(sig_values_text, cex=0.8)
  }
  dev.off()
}

# plot heatmaps for all the genes for all mags
met_gene_mag_datas <- list()
all_genes <- c(nit_genes, meth_genes)
for(i in 1:length(all_genes)){
  
  # prepare a common data frame for all MAGs for the marker gene
  temp_mark_gene_data <- data.frame(matrix(nrow = length(bin_mt_datas), ncol = ncol(all_mag_marker_gene_data)))
  rownames(temp_mark_gene_data) <- names(bin_mt_datas)
  colnames(temp_mark_gene_data) <- colnames(all_mag_marker_gene_data)
  
  for(j in 1:length(bin_mt_datas)){
    
    # get MAG data
    temp <- bin_mt_datas[[j]]
    temp <- temp$tpm$metmarkdb
    
    # parse sample names and organize similarly to other datas
    colnames(temp) <- paste("P",unlist(lapply(colnames(temp), function(x) strsplit(x = x, split = "-")[[1]][1])), sep = "")
    temp <- temp[,colnames(all_mag_marker_gene_data)]
    
    # check if the gene was detected for the MAG
    if(all_genes[i]%in%rownames(temp)){
      temp_mark_gene_data[j,] <- as.numeric(temp[all_genes[i],])
    }else{
      temp_mark_gene_data[j,] <- rep(0, ncol(temp_mark_gene_data))
    }
  }
  
  met_gene_mag_datas[[i]] <- temp_mark_gene_data 
}
names(met_gene_mag_datas) <- all_genes

# load prepared taxonomy for plotting 
load(paste(project_root, "/downstream/manuscript_figures/figure6/modified_gtdb_tax_silva.RData", sep = ""))

# prepared row names similarly as before for the other heatmaps
m_num <- seq(1:nrow(gtdbtk_report))
r_names <- paste("MAG",m_num, "; ",gtdbtk_report$Phylum,"; ", gtdbtk_report$Order, "; ", gtdbtk_report$Genus, sep = "")

# annotations for plotting 
annotation_col <- meta_use
for(j in 1:ncol(annotation_col)){annotation_col[,j] <- as.character(annotation_col[,j])}

annotation_colors <- list(
  grazing = c("grazed"="tan1", "ungrazed"="forestgreen"),
  veg_clusters = c("t_ces"="darkolivegreen2", "c_ros"="purple", "c_cho"="deeppink"),
  treatment = c("-S"="dodgerblue2", "CTL"="gainsboro", "+S"="firebrick3")
)

# plot the marker gene expression mag heatmaps
for(i in 1:length(met_gene_mag_datas)){
  data_heat <- met_gene_mag_datas[[i]]
  data_heat <- data_heat[rownames(gtdbtk_report),]
  rownames(data_heat) <- r_names
  
  # get good range for colors
  r <- c(as.numeric(unlist(data_heat)))
  r <- c(quantile(r, 0.05, na.rm = T), quantile(r, 0.995))
  if(r[2]==0){r[2] <- max(as.numeric(unlist(data_heat)))}
  
  # define colors
  breakslist <- seq(from=r[1], to=r[2], length.out=11)
  cols_heat <- rev(brewer.pal(11,"RdBu"))
  
  # define filename
  file_n <- paste(names(met_gene_mag_datas)[i], "_mag_heatmap.pdf", sep = "")
  
  # filter away all zero rows - no need to plot MAGs with no expression
  data_heat <- data_heat[-which(rowSums(data_heat)==0),]
  
  # clustering is not often possible due to lot's of zeroes, don't do
  # plot heatmap into object
  ph <- pheatmap(mat = data_heat, na_col = "darkgrey", cluster_rows = T, clustering_distance_rows = "correlation", fontsize = 12, breaks = breakslist, color = cols_heat, cluster_cols = F, annotation_col = annotation_col, annotation_colors = annotation_colors, scale="none", cellwidth = 14, cellheight = 14, silent = T, fontsize_row = 14, fontsize_col = 14, filename = NA)
  
  # italicize names
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
  pdf(file_n, width = 16, height = 10)
  grid.newpage()
  grid.draw(gt)
  dev.off()
}
