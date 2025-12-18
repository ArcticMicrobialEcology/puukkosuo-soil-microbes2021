# a small script to plot the marginal means for the LMM models for the methane and nitrogen related environmental variables

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(lmerTest)
library(emmeans)

# define custom functions
plotMarginalMeans <- function(emm_temp, ylab, xlab, cols, main_tit){
  # convert to a data frame if it's not already (confint.emmGrid typically returns a data.frame)
  ci_data <- as.data.frame(emm_temp)
  
  # set up the plot margins (optional, but can help with labels)
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # create an empty plot to define the axes
  plot(
    x = 1:nrow(ci_data), # use numerical indices for x-axis positioning
    y = ci_data$emmean,
    ylim = range(c(ci_data$lower.CL, ci_data$upper.CL)), # set y-axis limits to encompass CIs
    type = "n", # don't plot any points yet
    xaxt = "n", # suppress default x-axis
    xlab = xlab,
    ylab = ylab,
    main = main_tit,
    cex.lab = 1.2, # adjust label size
    cex.axis = 1.1, # adjust axis tick label size
    cex.main = 1.3, # adjust title size
    bty = "n" # remove box
  )
  
  # Add the confidence intervals as vertical lines
  segments(
    x0 = 1:nrow(ci_data),
    y0 = ci_data$lower.CL,
    x1 = 1:nrow(ci_data),
    y1 = ci_data$upper.CL,
    col = "black", # line color
    lwd = 2 # line width
  )
  
  # add vertical lines for the means
  for(i in 1:nrow(ci_data)){
    abline(h = ci_data$emmean[i], lty=2, lwd=2, col="grey")
  }
  
  # add the points for the estimated marginal means
  points(
    x = 1:nrow(ci_data),
    y = ci_data$emmean,
    pch = 19, # solid circle
    col = cols, # point color
    cex = 2 # point size
  )
  
  # add custom x-axis labels
  axis(
    side = 1, # bottom axis
    at = 1:nrow(ci_data), # positions for the labels
    labels = ci_data$grazing, # the actual grazing factor levels
    las = 1 # make labels perpendicular to the axis (optional)
  )
}

# load parsed gas flux and pore water data
load(paste(project_root, "/metadata/environmental_data/Meth_fluxes_porewater.RData", sep = ""))

# load parsed gas flux and pore water data
load(paste(project_root, "/downstream/manuscript_figures/figure5/All_meth_nit_models.RData", sep = ""))

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# estimate the marginal means for grazing for methane fluxes

# first date
emm_meth1 <- emmeans(mod_meth_flux1, ~ grazing)

# second date
emm_meth2 <- emmeans(mod_meth_flux2, ~ grazing)

# nitrate+nitrite
emm_nitrate_nitrite <- emmeans(mod_nitrate_nitrite, ~ grazing)

# total nitrogen
emm_total_nit <- emmeans(mod_total_nit, ~ grazing)

# plot the marginal means plots
# estimated marginal means of methane/nitrite+nitrate/total nitrogen by the exlusion treatment / grazing

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

{
  pdf(file = "marginal_means_methane_nitrogen.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plotMarginalMeans(emm_temp = emm_meth1, ylab = "Estimated CH4 FLUX [mg m-2 h-1]", xlab = "", cols = c("forestgreen", "tan1"), main_tit = "12_9_2021")
  plotMarginalMeans(emm_temp = emm_meth2, ylab = "Estimated CH4 FLUX [mg m-2 h-1]", xlab = "", cols = c("forestgreen", "tan1"), main_tit = "20_9_2021")
  plotMarginalMeans(emm_temp = emm_nitrate_nitrite, ylab = "Estimated NO3+NO2 µg/l", xlab = "", cols = c("forestgreen", "tan1"), main_tit = "NO3+NO2 µg/l")
  plotMarginalMeans(emm_temp = emm_total_nit, ylab = "Estimated TOTN µg/l", xlab = "", cols = c("forestgreen", "tan1"), main_tit = "TOTN µg/l")
  dev.off()
  
}
