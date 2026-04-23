# load libraries
library(readxl)

# define custom functions

# a small functio to get mean and sd
mean_sd <- function(x) {
  x <- as.numeric(x)
  c(
    mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE)
  )
}

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# read in pore-water datas
pore_water <- data.frame(read_excel(paste(project_root, "/metadata/environmental_data/Huokosvesi_Puukkosuo_kesä2021.xlsx", sep = ""), sheet = 1), stringsAsFactors = F, check.names = F)

# seperate months
pore_water_may <- pore_water[pore_water$Kuukausi=="toukokuu",]
pore_water_june <- pore_water[pore_water$Kuukausi=="kesäkuu",]
pore_water_july <- pore_water[pore_water$Kuukausi=="heinäkuu",]
pore_water_august <- pore_water[pore_water$Kuukausi=="elokuu",]
pore_water_september <- pore_water[pore_water$Kuukausi=="syyskuu",]

# read-in prepared metadata
load(paste(project_root,"/metadata/Study_Metadata.RData", sep = ""))
rownames(metadata) <- metadata$Plot

# read-in metadata
setwd("C:/Users/03261091/OneDrive - Luonnonvarakeskus/Oulanka_MG_MT_project_2009164/Metadata")
sample_table <- data.frame(read_excel("Oulanka-ACAP-study site.xlsx"), stringsAsFactors = F, check.names = F)
rownames(sample_table) <- sample_table$Plot

# organize according to metadata
rownames(pore_water_may) <- pore_water_may$RuutuID
rownames(pore_water_june) <- pore_water_june$RuutuID
rownames(pore_water_july) <- pore_water_july$RuutuID
rownames(pore_water_august) <- pore_water_august$RuutuID
rownames(pore_water_september) <- pore_water_september$RuutuID

pore_water_may <- pore_water_may[rownames(metadata),]
pore_water_june <- pore_water_june[rownames(metadata),]
pore_water_july <- pore_water_july[rownames(metadata),]
pore_water_august <- pore_water_august[rownames(metadata),]
pore_water_september <- pore_water_september[rownames(metadata),]

# add metadata to pore water tables
pore_water_may <- cbind(pore_water_may, metadata)
pore_water_june <- cbind(pore_water_june, metadata)
pore_water_july <- cbind(pore_water_july, metadata)
pore_water_august <- cbind(pore_water_august, metadata)
pore_water_september <- cbind(pore_water_september, metadata)

# two measurement values for plot number 69 in for august
# take the mean of these and use that, would something else be better?
double_input <- pore_water_august["69",5:11]
# first, replace , with . for easier parsing
double_input <- gsub(",",".", double_input)

# take the mean for each value
double_input <- unlist(lapply(strsplit(x = double_input, split = "/"), function(x) mean(as.numeric(x))))

# put averaged values back in to the august pore water table
pore_water_august["69",5:11] <- double_input

# now, put all months back into one table
pore_water <- rbind(pore_water_may, pore_water_june, pore_water_july, pore_water_august,pore_water_september)

# set the relevant variables as numeric
colnames(pore_water)
for(i in 5:11){pore_water[,i] <- as.numeric(pore_water[,i])}

# set also relevant variables as factors, missing
pore_water$Grazing <- factor(pore_water$Grazing, levels = c("grazed", "ungrazed"))
pore_water$Treatment <- factor(pore_water$Treatment, levels = c("CTL", "-S", "+S"))
pore_water$Microsite <- factor(pore_water$Microsite, levels = c("B.nana", "Lawn", "Wet"))
pore_water$Veg_cluster <- factor(pore_water$Veg_cluster, levels = c("T.Ces", "C.Cho", "C.Ros"))

# separate by the exclusion treatment
pore_water_grazed <- pore_water[pore_water$Grazing=="grazed",]
pore_water_ungrazed <- pore_water[pore_water$Grazing=="ungrazed",]

# get months, already in the correct order
months <- unique(pore_water_grazed$Kuukausi)

# in finnish, change to english in plotting
months_plot <- c("May", "June", "July", "August", "September")

# plot everything into file like this
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

# plot mean +- sd
{
  pdf(file = "pore_water_chemistry_summer_2021.pdf", width = 12, height = 8.3, onefile = T)
  par(mfrow=c(2,4))
  par(mar = c(8.1, 4.1, 4.1, 2.1))
  for(i in seq(from=5, to=11)){
    
    # get mean and sd for outside plots for the pore water variable
    grazed_stats <-data.frame(t(sapply(months, function(m) {
      mean_sd(pore_water_grazed[pore_water_grazed$Kuukausi == m, i])
    })), stringsAsFactors = FALSE, check.names = FALSE)
    
    # inside plots
    ungrazed_stats <- data.frame(t(sapply(months, function(m) {
      mean_sd(pore_water_ungrazed[pore_water_ungrazed$Kuukausi == m, i])
    })), stringsAsFactors = FALSE, check.names = FALSE)
    
    # calculate p-value with wilcoxon text
    p_vals <- sapply(months, function(m){
      wilcox.test(x = as.numeric(pore_water_grazed[pore_water_grazed$Kuukausi == m, i]), 
                  y=as.numeric(pore_water_ungrazed[pore_water_ungrazed$Kuukausi == m, i]), 
                  alternative = "two.sided")$p.value
    })
    
    # calculate a range for the plot, use standard deviation to show variability for each date.
    y_lim <- range(
      grazed_stats$mean - grazed_stats$sd,
      grazed_stats$mean + grazed_stats$sd,
      ungrazed_stats$mean - ungrazed_stats$sd,
      ungrazed_stats$mean + ungrazed_stats$sd,
      na.rm = TRUE
    )
    
    # plot
    x <- seq_along(months)
    
    # outside mean line
    plot(
      x, grazed_stats[, "mean"],
      type = "l", lwd = 3, col = "tan1",
      ylim = y_lim,
      xaxt = "n",
      xlab = "", ylab = colnames(pore_water_grazed)[i],
      main = "",
      cex = 1.5, cex.axis = 1.4, cex.lab = 1.5
    )
    
    # SD ribbon (outside)
    polygon(
      c(x, rev(x)),
      c(grazed_stats[, "mean"] - grazed_stats[, "sd"],
        rev(grazed_stats[, "mean"] + grazed_stats[, "sd"])),
      col = adjustcolor("tan1", alpha.f = 0.3),
      border = NA
    )
    
    # inside mean line
    lines(x, ungrazed_stats[, "mean"], lwd = 3, col = "forestgreen")
    
    # SD ribbon (inside)
    polygon(
      c(x, rev(x)),
      c(ungrazed_stats[, "mean"] - ungrazed_stats[, "sd"],
        rev(ungrazed_stats[, "mean"] + ungrazed_stats[, "sd"])),
      col = adjustcolor("forestgreen", alpha.f = 0.3),
      border = NA
    )
    
    # add months as axis
    axis(1, at = x, labels = months_plot, las = 2, cex = 1.5, cex.axis = 1.4, cex.lab = 1.5)
    
    # add annotation
    legend(
      "topright",
      legend = c("Outside", "Inside"),
      col = c("tan1", "forestgreen"),
      lwd = 3, bty = "n"
    )
    
    # add p-vals
    p_vals <- as.character(round(p_vals, 3))
    if(any(p_vals=="0")){p_vals[which(p_vals=="0")] <- "<0.001"}
    axis(3, at = x, labels = p_vals, las = 2, cex = 1.5, cex.axis = 1.4, cex.lab = 1.5)
  }
  dev.off()
}
