# a small script to plot all the pore water variables

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load parsed gas flux and pore water data
load(paste(project_root, "/metadata/environmental_data/Meth_fluxes_porewater.RData", sep = ""))

# take only the related september pore water variables
pore_water_september <- pore_water[which(pore_water$Kuukausi=="syyskuu"),]

# organize the pore water data similarly to the gas flux data, 
# according to the sample plot numberin in sample metadata
rownames(pore_water_september) <- pore_water_september$RuutuID
pore_water_september <- pore_water_september[rownames(metadata),]

# set the relevant variables as numeric
colnames(pore_water_september)
for(i in 5:11){pore_water_september[,i] <- as.numeric(pore_water_september[,i])}

# explore and plot the variables
# define the colors to be used for the grazing treatment
col_ungrazed <- "forestgreen"
col_grazed <- "tan1"

# set the levels for plotting differently
metadata$graz <- factor(metadata$Grazing, levels = c("ungrazed", "grazed"))

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

pdf(file = "pore_water_variables_september_exclusion_grazing.pdf", width = 14, height = 12, onefile = T)
par(mfrow=c(3,3))
par(mar=c(10.1, 4.1, 4.1, 2.1))
ylab <- "Value - unit in title"
for(i in 5:11){
  vals <- pore_water_september[,i]
  g <- as.factor(metadata$graz)
  boxplot(vals ~ g, ylab = ylab, xlab = "", las = 2, main = colnames(pore_water_september)[i], col = c(col_ungrazed, col_grazed), names=c("Ungrazed downslope", "Grazed upslope"),outline = F)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], cex=1.5, pch=19)
  }
  sig_val <- wilcox.test(x=vals[which(g=="grazed")], y=vals[which(g=="ungrazed")], alternative = "two.sided")
  mtext(paste("Wilcoxon test significance value:", round(sig_val$p.value, 3)))
}
dev.off()

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures