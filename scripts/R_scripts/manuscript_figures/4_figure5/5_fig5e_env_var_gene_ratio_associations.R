# a small script to plot the assocations of nitrogen and methane related environmental variabels to the related
# marker genes and ratios in MG and MT

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load parsed gas flux and pore water data
load(paste(project_root, "/metadata/environmental_data/Meth_fluxes_porewater.RData", sep = ""))

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

# define the specific genes investigated
nitrogen_genes <- c("NirK","NirS", "NorB", "NosZ")
methane_genes <- c("PmoA", "McrA")

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure5", sep = ""))

# start analysing and plotting
keep_items <- ls()
keep_items <- c(keep_items, "keep_items")
for(dat in 1:length(datas_list)){
  print(paste("Processing data:", names(datas_list)[dat]))
  
  # pick the relevant data
  data <- datas_list[[dat]]
  
  # nitrogen genes
  nit_genes <- nitrogen_genes[which(nitrogen_genes%in%rownames(data))]
  nit_data <- data[nit_genes,]
  
  # nirK+nirS/nosZ)
  nit_ratio2 <- (as.numeric(nit_data["NirK",]) +  as.numeric(nit_data["NirS",])) / as.numeric(nit_data["NosZ",])
  
  # methane genes
  meth_genes <- methane_genes[which(methane_genes%in%rownames(data))]
  meth_data <- data[meth_genes,]
  
  # ratio
  # pmoA / mcrA
  meth_ratio <- as.numeric(meth_data["PmoA",]) / as.numeric(meth_data["McrA",])
  
  # plot
  cols_use <- rep("tan1", nrow(meta_use))
  cols_use[meta_use$grazing=="ungrazed"] <- "forestgreen"
  
  # methane associated
  # plot with modified gene names - same as rownames for meth_data
  plot_names <- c("pmoA","mcrA")
  {
    pdf(file = paste(names(datas_list)[dat],"_methane_flux_genes_ratios_linear_associations.pdf"), width = 11.7, height = 7, onefile = T)
    par(mfrow=c(2,4))
    
    # first methane fluxes and ratio
    # ratio
    if(length(meth_ratio)>1){
      
      # log transform ratios also here
      ratio_vals <- log(meth_ratio)
      
      # log transform methane flux values
      flux_vals <- log(as.numeric(methane_fluxes[,1]))
      
      # use a simple linear model here to calculate associations (equivalet to Pearson correlation)
      mod_temp <- lm(flux_vals~ratio_vals)
      sum_temp <- summary(mod_temp)
      
      main_text <- "pmoA / mcrA"
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      plot(ratio_vals, flux_vals, pch=19, col=cols_use, ylab = "Ln-transformed CH4 FLUX [mg m-2 h-1]", xlab = " Ln-transformed PmoA / McrA", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5)
      abline(mod_temp, col="red", lwd=2)
      mtext(paste("P-value:", round(sum_temp$coefficients["ratio_vals",4],3)))
    }
    
    # genes
    for(i in 1:nrow(meth_data)){
      # log2 transform gene tpm
      ratio_vals <- log2(as.numeric(meth_data[i,])+1)
      
      # use a simple linear model here to calculate associations (equivalet to Pearson correlation)
      mod_temp <- lm(flux_vals~ratio_vals)
      sum_temp <- summary(mod_temp)
      
      main_text <- plot_names[i]
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      plot(ratio_vals, flux_vals, pch=19, col=cols_use, ylab = "Ln-transformed CH4 FLUX [mg m-2 h-1]", xlab = "Log2-transformed copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5)
      abline(mod_temp, col="red", lwd=2)
      mtext(paste("P-value:", round(sum_temp$coefficients["ratio_vals",4],3)))
    }
    
    # empty plot
    plot(1,1)
    
    # second methane flux, different data
    # ratio
    if(length(meth_ratio)>1){
      
      # log transform ratios also here
      ratio_vals <- log(meth_ratio)
      flux_vals <- log(as.numeric(methane_fluxes[,2]))
      
      # use a simple linear model here to calculate associations (equivalet to Pearson correlation)
      mod_temp <- lm(flux_vals~ratio_vals)
      sum_temp <- summary(mod_temp)
      
      main_text <- "pmoA / mcrA"
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      plot(ratio_vals, flux_vals, pch=19, col=cols_use, ylab = "Ln-transformed CH4 FLUX [mg m-2 h-1]", xlab = " Ln-transformed PmoA / McrA", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5)
      abline(mod_temp, col="red", lwd=2)
      mtext(paste("P-value:", round(sum_temp$coefficients["ratio_vals",4],3)))
    }
    
    # genes
    for(i in 1:nrow(meth_data)){
      
      # log2 transform gene values
      ratio_vals <- log2(as.numeric(meth_data[i,])+1)
      
      # use a simple linear model here to calculate associations (equivalet to Pearson correlation)
      mod_temp <- lm(flux_vals~ratio_vals)
      sum_temp <- summary(mod_temp)
      
      main_text <- plot_names[i]
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      plot(ratio_vals, flux_vals, pch=19, col=cols_use, ylab = "Ln-transformed CH4 FLUX [mg m-2 h-1]", xlab = "Log2-transformed copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5)
      abline(mod_temp, col="red", lwd=2)
      mtext(paste("P-value:", round(sum_temp$coefficients["ratio_vals",4],3)))
    }
    dev.off()
  }
  
  # nitrogen associated
  nit_data <- nit_data[c("NirS", "NorB", "NosZ"),]
  
  # plot with modified gene names - same as rownames for nit_data
  plot_names <- c("nirS","norB", "nosZ")
  
  {
    pdf(file = paste(names(datas_list)[dat], "_nitrogen_pore_water_genes_ratios_linear_associations.pdf"), width = 11.7, height = 7, onefile = T)
    par(mfrow=c(2,4))
    
    # nitrate + nitrite and ratio
    pore_water_vals <- log(as.numeric(pore_water_nitrogen$`NO3+NO2 µg/l`))
    if(length(nit_ratio2)>1){
      
      # log transform ratios also here
      ratio_vals <- log(nit_ratio2)
      
      # use a simple linear model here to calculate associations (equivalet to Pearson correlation)
      mod_temp <- lm(pore_water_vals~ratio_vals)
      sum_temp <- summary(mod_temp)
      
      main_text <- "(nirK + nirS) / nosZ"
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      plot(ratio_vals, pore_water_vals, pch=19, col=cols_use, ylab = "Ln-transformed NO3+NO2 µg/l", xlab = " Ln-transformed (NirK+NirS) / Nosz", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5)
      abline(mod_temp, col="red", lwd=2)
      mtext(paste("P-value:", round(sum_temp$coefficients["ratio_vals",4],3)))
    }
    
    # genes
    for(i in 1:nrow(nit_data)){
      
      # log2 transform gene values
      ratio_vals <- log2(as.numeric(nit_data[i,])+1)
      
      # use a simple linear model here to calculate associations (equivalet to Pearson correlation)
      mod_temp <- lm(pore_water_vals~ratio_vals)
      sum_temp <- summary(mod_temp)
      
      main_text <- plot_names[i]
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      plot(ratio_vals, pore_water_vals, pch=19, col=cols_use, ylab = "Ln-transformed NO3+NO2 µg/l", xlab = "Log2-transformed copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5)
      abline(mod_temp, col="red", lwd=2)
      mtext(paste("P-value:", round(sum_temp$coefficients["ratio_vals",4],3)))
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

