# a small script to calculate and plot the correlations between nitrogen and methane related marker gene ratios and taxa abundance and expression
# at MG and MT at the order level

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(phyloseq)
library(microbiome)

# calculate correlations between marker genes of interest and taxonomic data

# load the marker gene datas

# put datas into a list
datas_list <- list()

# metagenomics
load(paste(project_root, "/metagenomics/metmarkdb_diamond/Matrices_For_Downstream.RData", sep = ""))
datas_list[[1]] <- metmark_tpm_data

# metatranscriptomics
load(paste(project_root, "/metatranscriptomics/metmarkdb_diamond/Matrices_For_Downstream.RData", sep = ""))
datas_list[[2]] <- metmark_tpm_data

names(datas_list) <- c("metmark_mg", "metmark_mt")

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# load order level taxonomy data
# metagenomics
load(paste(project_root, "/metagenomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = ""))

# here we have compositional data already aggregated to the order level
psq_order_mg <- psq_order

# metatranscriptomics 
load(paste(project_root, "/metatranscriptomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = ""))

psq_order_mt <- psq_order

psq_list <- list(psq_order_mg, psq_order_mt)
names(psq_list) <- c("order_mg", "order_mt")

# clr transform and get compostional / relative NTU/OTU data
datas_clr <- list()
datas_rel <- list()

# define nitrogen and methane related genes
nitrogen_genes <- c("NirS", "NirK", "NorB", "NosZ")
methane_genes <- c("PmoA", "McrA")

for(i in 1:length(psq_list)){
  temp_psq <- psq_list[[i]]
  data_rel_temp <- data.frame(otu_table(temp_psq), stringsAsFactors = F, check.names = F)
  data_rel_temp <- data_rel_temp * 100
  datas_rel[[i]] <- data_rel_temp
  
  temp_psq <-  microbiome::transform(x = temp_psq, transform = "clr")
  
  datas_clr[[i]] <- data.frame(otu_table(temp_psq), stringsAsFactors = F, check.names = F)
}
names(datas_clr) <- names(datas_rel) <- names(psq_list)

# note: samples always organized in similar order in all datas
# calculate and plot correlations

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

# start analysing and plotting
keep_items <- ls()
keep_items <- c(keep_items, "keep_items")

for(dat in 1:length(datas_list)){
  
  print(paste("Processing data:", names(datas_list)[dat]))
  
  # pick the relevant data
  data <- datas_list[[dat]]
  data_clr <- datas_clr[[dat]]
  data_rel <- datas_rel[[dat]]
  
  # get otu mean abundances
  otu_mean_abundances <- rowMeans(data_rel)
  otu_mean_abundances <- sort(otu_mean_abundances, decreasing = T)
  
  # nirK+nirS/nosZ)
  nit_ratio2 <- (as.numeric(data["NirK",]) +  as.numeric(data["NirS",])) / as.numeric(data["NosZ",])
  
  # natural logarithm tranformation for the ratios here
  nit_ratio2 <- log(nit_ratio2)
  
  # pmoA / mcrA
  meth_ratio <- as.numeric(data["PmoA",]) / as.numeric(data["McrA",])
  
  # natural logarithm tranformation
  meth_ratio <- log(meth_ratio)
  
  # calculate correlations
  # nitrogen
  nit_cors <- data.frame(matrix(nrow = nrow(data_clr), ncol = 2), stringsAsFactors = F, check.names = F)
  for(ro in 1:nrow(data_clr)){
    temp_tax <- as.numeric(data_clr[ro,])
    cor_temp <- cor.test(nit_ratio2, temp_tax, method = "pearson", use = "pairwise.complete.obs")
    nit_cors[ro,1] <- cor_temp$estimate
    nit_cors[ro,2] <- cor_temp$p.value
  }
  rownames(nit_cors) <- rownames(data_clr)
  colnames(nit_cors) <- c("cor", "pvalue")
  
  # adjust p-values for multiple hypothesis testing
  nit_cors$fdr <- p.adjust(p = as.numeric(nit_cors$pvalue),method = "fdr")
  
  # methane
  meth_cors <- data.frame(matrix(nrow = nrow(data_clr), ncol = 2), stringsAsFactors = F, check.names = F)
  for(ro in 1:nrow(data_clr)){
    temp_tax <- as.numeric(data_clr[ro,])
    cor_temp <- cor.test(meth_ratio, temp_tax, method = "pearson", use = "pairwise.complete.obs")
    meth_cors[ro,1] <- cor_temp$estimate
    meth_cors[ro,2] <- cor_temp$p.value
  }
  rownames(meth_cors) <- rownames(data_clr)
  colnames(meth_cors) <- c("cor", "pvalue")
  meth_cors$fdr <- p.adjust(p = as.numeric(meth_cors$pvalue),method = "fdr")
  
  # plot some examples of best correlations
  # nitrogen
  
  # get only the significant correlations
  sig_cors <- nit_cors[which(nit_cors$fdr<=0.05),]
  
  # order significant correlations based on relative abundance of the taxa
  sig_abund <- otu_mean_abundances[rownames(sig_cors)]
  sig_abund <- sort(sig_abund, decreasing = T)
  sig_cors <- sig_cors[names(sig_abund),]
  
  # only three signficant correlations for MG
  # plot at most top 4 correlations for nitrogen related marker gene ratios
  if(nrow(sig_cors)>3){
    sig_cors <- sig_cors[1:4,]
  }
  
  # define plotting colors
  cols_use <- rep("tan1", nrow(meta_use))
  cols_use[which(meta_use$grazing=="ungrazed")] <- "forestgreen"
  
  {
    pdf(file = paste(names(datas_list)[dat], "_nitrogen_ratio_vs_top_taxa_examples.pdf", sep = ""), width = 11.7, height = 7, onefile = T)
    
    par(mfrow=c(2,4))
    for(ro in 1:nrow(sig_cors)){
      
      temp_tax <- as.numeric(data_clr[rownames(sig_cors)[ro],])
      
      main_text <- paste("(nirK + nirS) / nosZ vs.", rownames(sig_cors)[ro])
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      plot(temp_tax, nit_ratio2, pch=19, col=cols_use, ylab = "Ln-transformed (nirK + nirS) / nosZ", xlab = "Centered log ratio tranformed abundance", main = main_text, cex.main = 1.2, cex.axis = 1.4, cex.lab = 1)
      mtext(paste("Cor:", round(sig_cors$cor[ro],3), ", Mean abundance:", round(sig_abund[ro],2), "%"), cex = 0.8, line = 0.7)
      mtext(paste("P-value:", round(sig_cors$pvalue[ro],3), "False discovery rate:", round(sig_cors$fdr[ro],3)), cex = 0.8, line = -0.1)
    }
    dev.off()
  }
  
  # methane
  # get only the significant correlations
  sig_cors <- meth_cors[which(meth_cors$fdr<=0.05),]
  
  # order significant correlations based on relative abundance of the taxa
  sig_abund <- otu_mean_abundances[rownames(sig_cors)]
  sig_abund <- sort(sig_abund, decreasing = T)
  sig_cors <- sig_cors[names(sig_abund),]
  # too many significant correlations still to plot
  
  # plot the best correlations only
  sig_cors <- sig_cors[order(as.numeric(sig_cors$pvalue)),]
  
  # plot only the two best correlations as an example
  sig_cors <- sig_cors[1:2,]
  
  {
    pdf(file = paste(names(datas_list)[dat], "_methane_ratio_vs_top_taxa_examples.pdf", sep = ""), width = 11.7, height = 7, onefile = T)
    
    par(mfrow=c(2,4))
    for(ro in 1:nrow(sig_cors)){
      
      temp_tax <- as.numeric(data_clr[rownames(sig_cors)[ro],])
      
      main_text <- paste("pmoA / mcrA vs. ", rownames(sig_cors)[ro])
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      plot(temp_tax, meth_ratio, pch=19, col=cols_use, ylab = "Ln-transformed pmoA / mcrA", xlab = "Centered log ratio tranformed abundance", main = main_text, cex.main = 1.2, cex.axis = 1.4, cex.lab = 1)
      mtext(paste("Cor:", round(sig_cors$cor[ro],3), ", Mean abundance:", round(sig_abund[ro],2), "%"), cex = 0.8, line = 0.7)
      mtext(paste("P-value:", round(sig_cors$pvalue[ro],3), "False discovery rate:", round(sig_cors$fdr[ro],3)), cex = 0.8, line = -0.1)
    }
    dev.off()
  }
  
  # clean
  del_items <- ls()
  del_items <- del_items[-which(del_items%in%keep_items)]
  rm(list = del_items)
  
  print("Data processed.")
  print("***************************************************************")
  
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures

