# a small script to plot the soil temperature and snow depth data

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(readxl)
library(lmerTest)

# change to environmental metadata directory
setwd(paste(project_root, "/metadata/environmental_data", sep = ""))

# read-in the the snow depth data
snow_depth_march <- data.frame(read_excel(path = "Temperature_Snow_depth.xlsx", sheet = 3), stringsAsFactors = F, check.names = F)

# define the colors to be used for the grazing treatment
col_ungrazed <- "forestgreen"
col_grazed <- "tan1"

# plot with respect to grazing
vals <- as.numeric(snow_depth_march$Snow_cm)
g <- factor(snow_depth_march$Grazing, levels = c("Ungrazed", "Grazed"))
labs <- c("Ungrazed downslope", "Grazed upslope")

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

# plot snow depth with respect to grazing
pdf(file = "snow_depth_grazing_march_2021.pdf", width = 8, height = 8)
par(mar = c(11.1, 4.1, 4.1, 2.1))
boxplot(vals ~ g, ylab = "Snow depth (cm)", xlab = "", names = labs, las = 2, outline = F, main = "Snow depth March 2021", col = c(col_ungrazed, col_grazed))
for(j in 1:length(unique(g))){
  points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], cex=1.5, pch=19)
}
sig_val <- wilcox.test(x=vals[which(g=="Grazed")], y=vals[which(g=="Ungrazed")], alternative = "two.sided")
mtext(paste("Wilcoxon test significance value:", round(sig_val$p.value, 3)))
dev.off()

# plot with respect to snow treatments
g <- factor(snow_depth_march$Treatment, levels = c("Ctl", "Add", "Rem"))
labs <- c("AMB", "+S", "-S")

pdf(file = "snow_depth_snow_treatment_march_2021.pdf", width = 8, height = 8)
par(mar = c(6.1, 4.1, 4.1, 2.1))
boxplot(vals ~ g, ylab = "Snow depth (cm)", xlab = "", names = labs, las = 2, outline = F, main = "Snow depth March 2021", col = c("gainsboro", "firebrick3", "dodgerblue2"))
for(j in 1:length(unique(g))){
  points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], cex=1.5, pch=19)
}
sig_val <- kruskal.test(vals~g)
mtext(paste("Kruskal - Wallis test significance value:", round(sig_val$p.value, 3)))
dev.off()

# soil temperature
# change to environmental metadata directory
setwd(paste(project_root, "/metadata/environmental_data", sep = ""))

# read-in the the soil temp data
soil_temp <- data.frame(read_excel(path = "Temperature_Snow_depth.xlsx", sheet = 1), stringsAsFactors = F, check.names = F)

# different annotation used here for the exclusion / grazing treatment, but doesn't matter, everything will postedited in inkscape
soil_temp$Grazing <- factor(soil_temp$Grazing, levels = c("UGr_downslope", "Gr_upslope"))
soil_temp$Snow <- factor(soil_temp$Snow, levels = c("Ctl", "Add", "Rem"))

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

# plot month by month
all_months <- unique(soil_temp$Month)
labs <- c("Ungrazed downslope", "Grazed upslope")

pdf(file = "soil_temp_grazing.pdf", width = 11, height = 12, onefile = T)
par(mfrow=c(4,3))
ylab <- "Temperature (C)"
for(i in 1:length(all_months)){
  par(mar = c(11.1, 4.1, 4.1, 2.1))
  temp <- soil_temp[which(soil_temp$Month==all_months[i]),]
  vals <- as.numeric(temp$Temp)
  g <- factor(temp$Grazing, levels = c("UGr_downslope", "Gr_upslope"))
  boxplot(vals ~ g, ylab = ylab, xlab = "", las = 2, names = labs, main = all_months[i], col = c(col_ungrazed, col_grazed), outline = F)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], cex=1.5, pch=19)
  }
  sig_val <- wilcox.test(x=vals[which(g=="Gr_upslope")], y=vals[which(g=="UGr_downslope")], alternative = "two.sided")
  mtext(paste("Wilcoxon test significance value:", round(sig_val$p.value, 3)))
}
dev.off()
# soil temp not different at any month.

# snow treatment
labs <- c("AMB", "+S", "-S")
pdf(file = "soil_temp_snow_treatment.pdf", width = 11, height = 12, onefile = T)
par(mfrow=c(4,3))
ylab <- "Temperature (C)"
for(i in 1:length(all_months)){
  temp <- soil_temp[which(soil_temp$Month==all_months[i]),]
  vals <- as.numeric(temp$Temp)
  g <- factor(temp$Snow, levels = c("Ctl", "Add", "Rem"))
  boxplot(vals ~ g, ylab = ylab, xlab = "", names = labs, las = 2, main = all_months[i], col = c("gainsboro", "firebrick3", "dodgerblue2"), outline = F)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], cex=1.5, pch=19)
  }
  sig_val <- kruskal.test(vals~g)
  mtext(paste("Kruskal - Wallis test significance value:", round(sig_val$p.value, 3)))
}
dev.off()

# print out session info
print("SessionInfo:")
sessionInfo()

