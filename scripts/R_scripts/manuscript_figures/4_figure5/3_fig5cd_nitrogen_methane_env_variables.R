# a small script to plot the nitrogen related pore water variables and methane fluxes and perform LMM analysis for them

# load libraries 
library(lmerTest)

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load parsed gas flux and pore water data
load(paste(project_root, "/metadata/environmental_data/Meth_fluxes_porewater.RData", sep = ""))

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# separate grazing into two classes based on vegetation
meta_use$grazing_vegetation <- factor(paste(meta_use$grazing, meta_use$veg_clusters, sep = "_"), levels=c("ungrazed_t_ces","grazed_t_ces", "ungrazed_c_cho", "grazed_c_cho", "grazed_c_ros"))

# model methane fluxes according to grazing
temp_data <- cbind(log(methane_fluxes[,1]), meta_use)
colnames(temp_data)[1] <- "methane"

# prepare both more complex model with interactions and a simpler model without interactions
set.seed(1)
mod_meth_flux1_simple <- lmer(methane ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)

set.seed(1)
mod_meth_flux1_complex <- lmer(methane ~ grazing * treatment + (1|veg_clusters),REML=FALSE, data = temp_data)

# compare the models
model_comparison <- anova(mod_meth_flux1_simple, mod_meth_flux1_complex)

# extract significance value
model_comparison_sig <- model_comparison$`Pr(>Chisq)`[2]

# save the the best model
if(model_comparison_sig<0.05){
  mod_meth_flux1 <- mod_meth_flux1_complex
  print("Best methane flux model is complex")
}else{
  mod_meth_flux1 <- mod_meth_flux1_simple
  print("Best methane flux model is simple")
}
# simple is the best model

# date 2
temp_data <- cbind(log(methane_fluxes[,2]), meta_use)
colnames(temp_data)[1] <- "methane"

# prepare both more complex model with interactions and a simpler model without interactions
set.seed(1)
mod_meth_flux2_simple <- lmer(methane ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)

set.seed(1)
mod_meth_flux2_complex <- lmer(methane ~ grazing * treatment + (1|veg_clusters),REML=FALSE, data = temp_data)

# compare the models
model_comparison <- anova(mod_meth_flux2_simple, mod_meth_flux2_complex)

# extract significance value
model_comparison_sig <- model_comparison$`Pr(>Chisq)`[2]

# save the the best model
if(model_comparison_sig<0.05){
  mod_meth_flux2 <- mod_meth_flux2_complex
  print("Best methane flux model is complex")
}else{
  mod_meth_flux2 <- mod_meth_flux2_simple
  print("Best methane flux model is simple")
}
# simple is the best model

# save the models into a list
meth_models <- list(mod_meth_flux1, mod_meth_flux2)

# plot methane fluxes with respect to grazing / exclusion treatment
# plot similarly to the marker genes

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure5", sep = ""))

# plot
{
  pdf(file = "env_variables_methane_fluxes_grazing.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  y_lab <- "CH4 FLUX [mg m-2 h-1]"
  for(i in 1:ncol(methane_fluxes)){
    vals <- methane_fluxes[,i]
    r_boxplot <- range(vals, na.rm = T)
    g <- meta_use$grazing
    boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, outline = F, ylab = y_lab, main = colnames(methane_fluxes)[i], cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, frame = FALSE)
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
    
    # significance value
    sig_val_temp <- summary(meth_models[[i]])
    sig_val_temp <- round(sig_val_temp$coefficients[2,5],3)
    sig_values_text <- paste("P-value:", sig_val_temp)
    mtext(sig_values_text, cex=0.8)
  }
  
  # methane fluxes with grazing divided accoring to vegetation
  for(i in 1:ncol(methane_fluxes)){
    vals <- methane_fluxes[,i]
    r_boxplot <- range(vals, na.rm = T)
    g <- meta_use$grazing_vegetation
    
    boxplot(vals~g, col=c("forestgreen","tan3","green3","tan2", "tan1"), xaxt="n", xlab = "", ylim = r_boxplot, outline = F, las=2, ylab = y_lab, main = colnames(methane_fluxes)[i], cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, frame = FALSE)
    axis(side = 1, at = c(1,2,3,4,5), labels = c("Ungrazed downslope T.Ces", "Grazed upslope T.Ces", "Ungrazed downslope C.Cho", "Grazed upslope C.cho", "Grazed upslope C.Ros"), las=2)
    
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
    
    # significance value
    sig_val_temp <- summary(meth_models[[i]])
    sig_val_temp <- round(sig_val_temp$coefficients[2,5],3)
    sig_values_text <- paste("P-value:", sig_val_temp)
    mtext(sig_values_text, cex=0.8)
  }
  dev.off() 
}

# pore water variables
# nitrate + nitrite
temp_data <- cbind(log(pore_water_nitrogen[,1]), meta_use)
colnames(temp_data)[1] <- "nitrate_nitrite"

# prepare both more complex model with interactions and a simpler model without interactions
set.seed(1)
mod_nitrate_nitrite_simple <- lmer(nitrate_nitrite ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)

set.seed(1)
mod_nitrate_nitrite_complex <- lmer(nitrate_nitrite ~ grazing * treatment + (1|veg_clusters),REML=FALSE, data = temp_data)

# compare the models
model_comparison <- anova(mod_nitrate_nitrite_simple, mod_nitrate_nitrite_complex)

# extract significance value
model_comparison_sig <- model_comparison$`Pr(>Chisq)`[2]

# save the the best model
if(model_comparison_sig<0.05){
  mod_nitrate_nitrite <- mod_nitrate_nitrite_complex
  print("Best nitrate+nitrite flux model is complex")
}else{
  mod_nitrate_nitrite <- mod_nitrate_nitrite_simple
  print("Best nitrate+nitrite flux model is simple")
}
# best model is simple

# total nitrogen
temp_data <- cbind(log(pore_water_nitrogen[,2]), meta_use)
colnames(temp_data)[1] <- "total_nit"

# prepare both more complex model with interactions and a simpler model without interactions
set.seed(1)
mod_total_nit_simple <- lmer(total_nit ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)

set.seed(1)
mod_total_nit_complex <- lmer(total_nit ~ grazing * treatment + (1|veg_clusters),REML=FALSE, data = temp_data)

# compare the models
model_comparison <- anova(mod_total_nit_simple, mod_total_nit_complex)

# extract significance value
model_comparison_sig <- model_comparison$`Pr(>Chisq)`[2]

# save the the best model
if(model_comparison_sig<0.05){
  mod_total_nit <- mod_total_nit_complex
  print("Best total nitrogen flux model is complex")
}else{
  mod_total_nit <- mod_total_nit_simple
  print("Best total nitrogen flux model is simple")
}
# best model is simple

# put the models in a list
nit_models <- list(mod_nitrate_nitrite, mod_total_nit)

{
  pdf(file = "env_variables_pore_water_nitrogen_grazing.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  for(i in 1:ncol(pore_water_nitrogen)){
    y_lab <- colnames(pore_water_nitrogen)[i]
    vals <- pore_water_nitrogen[,i]
    r_boxplot <- range(vals, na.rm = T)
    g <- meta_use$grazing
    boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, outline = F, ylab = y_lab, main = colnames(pore_water_nitrogen)[i], cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, frame = FALSE)
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
    
    # significance value
    sig_val_temp <- summary(nit_models[[i]])
    sig_val_temp <- round(sig_val_temp$coefficients[2,5],3)
    sig_values_text <- paste("P-value:", sig_val_temp)
    mtext(sig_values_text, cex=0.8)
  }
  
  
  # nitrogen with grazing divided into classes
  for(i in 1:ncol(pore_water_nitrogen)){
    y_lab <- colnames(pore_water_nitrogen)[i]
    vals <- pore_water_nitrogen[,i]
    r_boxplot <- range(vals, na.rm = T)
    g <- meta_use$grazing_vegetation
    
    boxplot(vals~g, col=c("forestgreen","tan3","green3","tan2", "tan1"), xaxt="n", xlab = "", ylim = r_boxplot, outline = F, las=2, ylab = y_lab, main = colnames(pore_water_nitrogen)[i], cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5, frame = FALSE)
    axis(side = 1, at = c(1,2,3,4,5), labels = c("Ungrazed downslope T.Ces", "Grazed upslope T.Ces", "Ungrazed downslope C.Cho", "Grazed upslope C.cho", "Grazed upslope C.Ros"), las=2)
    
    for(j in 1:length(unique(g))){
      points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
    }
    
    # significance value
    sig_val_temp <- summary(nit_models[[i]])
    sig_val_temp <- round(sig_val_temp$coefficients[2,5],3)
    sig_values_text <- paste("P-value:", sig_val_temp)
    mtext(sig_values_text, cex=0.8)
  }
  dev.off() 
}

# save models
ses_info <- sessionInfo()

save(mod_meth_flux1_simple, mod_meth_flux1_complex, mod_meth_flux1, mod_meth_flux2_simple, mod_meth_flux2_complex, mod_meth_flux2,
     mod_nitrate_nitrite_simple, mod_nitrate_nitrite_complex, mod_nitrate_nitrite, mod_total_nit_simple, mod_total_nit_complex,
     mod_total_nit, ses_info, file = "All_meth_nit_models.RData")

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures