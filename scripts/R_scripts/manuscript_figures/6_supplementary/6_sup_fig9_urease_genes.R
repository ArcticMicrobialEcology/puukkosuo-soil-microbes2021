# load libraries
library(lmerTest)

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load mg data
load(paste(project_root, "/metagenomics/kegg_diamond/Matrices_For_Downstream.RData", sep = ""))
kegg_mg <- ko_tpm_data

# load mt data
load(load(paste(project_root, "/metatranscriptomics/kegg_diamond/Matrices_For_Downstream.RData", sep = "")))
kegg_mt <- ko_tpm_data

# explore the urease subunits
urea_related <- c("K01428", "K01429", "K01430")
annots <- c("K01428; ureC; urease subunit alpha", "K01429; ureB; urease subunit beta", "K01430; ureA; urease subunit gamma" )

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# plot these
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))
# in MG data
{
  pdf(file = "urea_related_mg.pdf", width = 11.7, height = 8.3)
  par(mfrow=c(2,2))
  for(i in 1:length(urea_related)){
    if(!urea_related[i]%in%rownames(kegg_mg)){next}
    vals <- as.numeric(kegg_mg[urea_related[i],])
    
    g <- as.factor(meta_use$grazing)
    r_boxplot <- range(vals, na.rm = T)
    main_text <- annots[i]
    
    boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, ylab = "Copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
    
    # add significance values
    temp_data <- cbind(vals, meta_use)
    colnames(temp_data)[1] <- "vals"
    mod_urea <- lmer(vals ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
    sum_temp <- summary(mod_urea)
    mtext(paste("P-value:", round(sum_temp$coefficients[2,5],3)), cex=1.4)
  }
  
  urea_mg <- colSums(kegg_mg[urea_related,], na.rm = TRUE)
  vals <- urea_mg
  g <- as.factor(meta_use$grazing)
  r_boxplot <- range(vals, na.rm = T)
  main_text <- "Urease (K01428 + K01429 + K01430)"
  
  boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, ylab = "Copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  
  # add significance values
  temp_data <- cbind(vals, meta_use)
  colnames(temp_data)[1] <- "vals"
  mod_urea <- lmer(vals ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
  sum_temp <- summary(mod_urea)
  mtext(paste("P-value:", round(sum_temp$coefficients[2,5],3)), cex=1.4)
  dev.off()
}

# in MT data
{
  pdf(file = "urea_related_mt.pdf", width = 11.7, height = 8.3)
  par(mfrow=c(2,2))
  for(i in 1:length(urea_related)){
    if(!urea_related[i]%in%rownames(kegg_mt)){next}
    vals <- as.numeric(kegg_mt[urea_related[i],])
    
    g <- as.factor(meta_use$grazing)
    r_boxplot <- range(vals, na.rm = T)
    main_text <- annots[i]
    
    boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, ylab = "Copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
    
    # add significance values
    temp_data <- cbind(vals, meta_use)
    colnames(temp_data)[1] <- "vals"
    mod_urea <- lmer(vals ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
    sum_temp <- summary(mod_urea)
    mtext(paste("P-value:", round(sum_temp$coefficients[2,5],3)), cex=1.4)
  }
  
  urea_mt <- colSums(kegg_mt[urea_related,], na.rm = TRUE)
  vals <- urea_mt
  g <- as.factor(meta_use$grazing)
  r_boxplot <- range(vals, na.rm = T)
  main_text <- "Urease (K01428 + K01430)"
  
  boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, ylab = "Copies per million", main = main_text, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  
  # add significance values
  temp_data <- cbind(vals, meta_use)
  colnames(temp_data)[1] <- "vals"
  mod_urea <- lmer(vals ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
  sum_temp <- summary(mod_urea)
  mtext(paste("P-value:", round(sum_temp$coefficients[2,5],3)), cex=1.4)
  dev.off()
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures