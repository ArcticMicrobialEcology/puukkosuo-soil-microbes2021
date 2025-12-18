# a small script to plot the abundance and expression of the nitrogen and methane related metabolic marker genes,
# gene ratios and perform LMM analysis to the gene ratios. 

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(lmerTest)

# define nitrogen and methane related genes to be analyzed
nitrogen_genes <- c("NirS", "NirK", "NorB", "NosZ","NifH", "NapA", "NrfA", "NarG", "HzsA")
methane_genes <- c("PmoA", "McrA")

# put datas into a list
datas_list <- list()

# metabolic marker gene datas
# metagenomics
load(paste(project_root, "/downstream/contig_based/metagenomics/metmarkdb_diamond/Matrices_For_Downstream.RData", sep = ""))

# investigate the nitrogen and methane related genes
# nitrogen
data <- metmark_tpm_data
nit_genes <- nitrogen_genes[which(nitrogen_genes%in%rownames(data))]
nit_data <- data[nit_genes,]

# methane
meth_genes <- methane_genes[which(methane_genes%in%rownames(data))]
meth_data <- data[meth_genes,]

# combine into one
data <- rbind(nit_data, meth_data)

# add into list
datas_list[[1]] <- data

# metatranscriptomics
load(paste(project_root, "/downstream/contig_based/metatranscriptomics/metmarkdb_diamond/Matrices_For_Downstream.RData", sep = ""))

# nitrogen
data <- metmark_tpm_data
nit_genes <- nitrogen_genes[which(nitrogen_genes%in%rownames(data))]
nit_data <- data[nit_genes,]

# methane
meth_genes <- methane_genes[which(methane_genes%in%rownames(data))]
meth_data <- data[meth_genes,]

# combine into one
data <- rbind(nit_data, meth_data)

# add into list
datas_list[[2]] <- data

# name the datas
names(datas_list) <- c("metmark_mg", "metmark_mt")

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# separate grazing / exclusion treatment based on vegetation
meta_use$grazing_vegetation <- factor(paste(meta_use$grazing, meta_use$veg_clusters, sep = "_"), levels=c("ungrazed_t_ces","grazed_t_ces", "ungrazed_c_cho", "grazed_c_cho", "grazed_c_ros"))

# go throgh the MG and MT datas similarly
keep_items <- ls()
keep_items <- c(keep_items, "keep_items")
for(dat in 1:length(datas_list)){
  
  print(paste("Processing data:", names(datas_list)[dat]))
  
  # load the lmm results
  setwd(paste(project_root, "/downstream/contig_based/lmm_analysis", sep = ""))

  setwd(names(datas_list)[dat])
  load("lmm_analysis_done.RData")
  
  # pick the relevant data
  data <- datas_list[[dat]]
  
  # start plotting and testing
  # grazing exclusion treatment
  # make sure data is in correct order
  data <- data[,c(rownames(meta_use))]
  
  # nitrogen genes
  nit_genes <- nitrogen_genes[which(nitrogen_genes%in%rownames(data))]
  nit_data <- data[nit_genes,]
  temp_sig <- sig_values_lmm_simple[nit_genes,]
  
  # marker gene ratios
  # construct similar LMM models as for the genes and get significance values
  
  
  if(dat==2){
    # nirS is not found in contig MT data, ratios not defined
    nit_ratio1 <- nit_ratio2 <- mod_nit_ratio1_simple <- mod_nit_ratio1_complex <- model_nit_ratio1 <- mod_nit_ratio2_simple <-  mod_nit_ratio2_complex <- model_nit_ratio2 <- NA
    
  } else {
    # nrfA - (nirK+nirS)
    nit_ratio1 <- as.numeric(nit_data["NrfA",]) - (as.numeric(nit_data["NirK",]) +  as.numeric(nit_data["NirS",]))
    
    # cannot log transform, negative values, don't log transform
    # temp_data <- cbind(log(nit_ratio1), meta_use)
    temp_data <- cbind(nit_ratio1, meta_use)
    
    set.seed(1)
    mod_nit_ratio1_simple <- lmer(nit_ratio1 ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
    
    set.seed(1)
    mod_nit_ratio1_complex <- lmer(nit_ratio1 ~ grazing * treatment + (1|veg_clusters),REML=FALSE, data = temp_data)
    
    model_comparison <- anova(mod_nit_ratio1_simple, mod_nit_ratio1_complex)
    
    # extract significance value
    model_comparison_sig <- model_comparison$`Pr(>Chisq)`[2]
    
    #save the the best model
    if(model_comparison_sig<0.05){
      model_nit_ratio1 <- mod_nit_ratio1_complex
      print("Best nrfA - (nirK+nirS) model is complex")
    }else{
      model_nit_ratio1 <- mod_nit_ratio1_simple
      print("Best nrfA - (nirK+nirS) model is simple")
    }
    names(nit_ratio1) <- colnames(data)
    # best model for both MG and MT is simple, without interactions
    
    # (nirK+nirS)/nosZ
    nit_ratio2 <- (as.numeric(nit_data["NirK",]) +  as.numeric(nit_data["NirS",])) / as.numeric(nit_data["NosZ",])
    if(any(is.infinite(nit_ratio2))){nit_ratio2[which(is.infinite(nit_ratio2))]=NA}
    #temp_data <- cbind(log(nit_ratio2), meta_use)
    temp_data <- cbind(nit_ratio2, meta_use)
    
    set.seed(1)
    mod_nit_ratio2_simple <- lmer(nit_ratio2 ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
    
    set.seed(1)
    mod_nit_ratio2_complex <- lmer(nit_ratio2 ~ grazing * treatment + (1|veg_clusters),REML=FALSE, data = temp_data)
    
    model_comparison <- anova(mod_nit_ratio2_simple, mod_nit_ratio2_complex)
    
    # extract significance value
    model_comparison_sig <- model_comparison$`Pr(>Chisq)`[2]
    
    #save the the best model
    if(model_comparison_sig<0.05){
      model_nit_ratio2 <- mod_nit_ratio2_complex
      print("Best (nirK+nirS)/nosZ model is complex")
    }else{
      model_nit_ratio2 <- mod_nit_ratio2_simple
      print("Best (nirK+nirS)/nosZ model is simple")
    }
    names(nit_ratio2) <- colnames(data)
    # best model for both MG and MT is simple, without interactions
  }
  
  # plot with modified gene names - same as rownames for nit_data
  if(dat==1){
    plot_names <- c("nirS","nirK","norB","nosZ","nifH","napA","nrfA","narG")
  } else {
    plot_names <- c("nirK","norB","nosZ","nifH","napA","nrfA","narG")
  }
  
   # change directory into the plotting directory - needs to exist
  setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))
  
  # boxplot several in one page
  {
    pdf(file = paste(names(datas_list)[dat], "_grazing_tpm_values_boxplot_all_one_page_nitrogen_contig.pdf", sep = ""), width = 11.7, height = 8.3, onefile = T)
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
      sig_values_text <- paste("P-value:",round(sig_values_lmm_simple[rownames(nit_data)[ro],1],3))
      mtext(sig_values_text, cex=0.8)
    }
    
    # ratios
    if(length(nit_ratio1)>1){
      vals <- as.numeric(nit_ratio1)
      r_boxplot <- range(vals, na.rm = T)
      main_text <- "nrfA - (nirK+nirS)"
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      g <- as.factor(meta_use$grazing)
      boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, outline = F, ylab = "nrfA - (nirK+nirS)", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5,  frame = FALSE)
      
      for(j in 1:length(unique(g))){
        points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
      }
      
      # significance value
      sig_val_temp <- summary(model_nit_ratio1)
      sig_val_temp <- round(sig_val_temp$coefficients[2,5],3)
      sig_values_text <- paste("P-value:", sig_val_temp)
      mtext(sig_values_text, cex=0.8)
      
    }
    
    if(length(nit_ratio2)>1){
      vals <- as.numeric(nit_ratio2)
      r_boxplot <- range(vals, na.rm = T)
      main_text <- "(nirK+nirS) / nosZ"
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      g <- as.factor(meta_use$grazing)
      boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, outline = F, ylab = "(nirK+nirS) / nosZ", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5,  frame = FALSE)
      
      for(j in 1:length(unique(g))){
        points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
      }
      
      # significance value
      sig_val_temp <- summary(model_nit_ratio2)
      sig_val_temp <- round(sig_val_temp$coefficients[2,5],3)
      sig_values_text <- paste("P-value:", sig_val_temp)
      mtext(sig_values_text, cex=0.8)
    }
    dev.off()
  }
  
  
  # methane genes
  meth_genes <- methane_genes[which(methane_genes%in%rownames(data))]
  meth_data <- data[meth_genes,]
  temp_sig <- sig_values_lmm_simple[meth_genes,]
  
  # ratio
  # PmoA / McrA
  meth_ratio <- as.numeric(meth_data["PmoA",]) / as.numeric(meth_data["McrA",])
  if(any(is.infinite(meth_ratio))){meth_ratio[which(is.infinite(meth_ratio))]=NA}
  
  #temp_data <- cbind(log(meth_ratio), meta_use)
  temp_data <- cbind(meth_ratio, meta_use)
  
  set.seed(1)
  mod_meth_ratio_simple <- lmer(meth_ratio ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
  
  set.seed(1)
  mod_meth_ratio_complex <- lmer(meth_ratio ~ grazing * treatment + (1|veg_clusters),REML=FALSE, data = temp_data)
  
  model_comparison <- anova(mod_meth_ratio_simple, mod_meth_ratio_complex)
  
  # extract significance value
  model_comparison_sig <- model_comparison$`Pr(>Chisq)`[2]
  
  # save the the best model
  if(model_comparison_sig<0.05){
    model_meth_ratio <- mod_meth_ratio_complex
    print("Best pmoA / mcrA model is complex")
  }else{
    model_meth_ratio <- mod_meth_ratio_simple
    print("Best pmoA / mcrA model is simple")
  }
  names(meth_ratio) <- colnames(data)
  
  # simple is the best model for MG and MT
  
  # plot with modified gene names - same as rownames for nit_data
  plot_names <- c("pmoA","mcrA")
  # several in one page
  {
    pdf(file = paste(names(datas_list)[dat],"_grazing_tpm_values_boxplot_all_one_page_methane_contig.pdf"), width = 11.7, height = 8.3, onefile = T)
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    par(mfrow=c(2,5))
    for(ro in 1:nrow(meth_data)){
      vals <- as.numeric(meth_data[ro,])
      r_boxplot <- range(vals, na.rm = T)
      main_text <- plot_names[ro]
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      g <- as.factor(meta_use$grazing)
      boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, outline = F, ylab = "Copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5,  frame = FALSE)
      
      for(j in 1:length(unique(g))){
        points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
      }
      
      # add sig values
      sig_values_text <- paste("P-value:",round(sig_values_lmm_simple[rownames(meth_data)[ro],1],3))
      mtext(sig_values_text, cex=0.8)
    }
    
    # ratios
    if(length(meth_ratio)>1){
      vals <- as.numeric(meth_ratio)
      r_boxplot <- range(vals, na.rm = T)
      main_text <- "pmoA / mcrA"
      # italize
      main_text <- bquote(italic(.(main_text)))
      
      g <- as.factor(meta_use$grazing)
      boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, outline = F, ylab = "pmoA / mcrA", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, frame = FALSE)
      
      for(j in 1:length(unique(g))){
        points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
      }
      
      # significance value
      sig_val_temp <- summary(model_meth_ratio)
      sig_val_temp <- round(sig_val_temp$coefficients[2,5],3)
      sig_values_text <- paste("P-value:", sig_val_temp)
      mtext(sig_values_text, cex=0.8)
    }
    dev.off()
  }
  
  # save the lmm models for ratios
  ses_info <- sessionInfo()
  save(mod_nit_ratio1_simple, mod_nit_ratio1_complex, model_nit_ratio1, 
       mod_nit_ratio2_simple, mod_nit_ratio2_complex, model_nit_ratio2,
       mod_meth_ratio_simple, mod_meth_ratio_complex, model_meth_ratio, ses_info, file = paste(names(datas_list)[dat], "_LMM_models_ratios_contig.RData", sep = ""))
  
  # clean
  del_items <- ls()
  del_items <- del_items[-which(del_items%in%keep_items)]
  rm(list = del_items)
  
  print("Data processed.")
  print("***************************************************************")
  
  setwd("../")
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures

