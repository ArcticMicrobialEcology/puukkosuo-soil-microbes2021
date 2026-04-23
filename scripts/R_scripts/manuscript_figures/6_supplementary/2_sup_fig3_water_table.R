# a small script to plot water table depth for the summer 2021

# load libraries
library(readxl)

# define custom functions
# a function to calculate means, standard deviations and standard errors
col_stats <- function(df) {
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  
  mean <- colMeans(mat, na.rm = TRUE)
  sd   <- apply(mat, 2, sd, na.rm = TRUE)
  n    <- apply(mat, 2, function(x) sum(!is.na(x)))
  se   <- sd / sqrt(n)
  
  list(mean = mean, sd = sd, se = se)
}

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# read in water-level data. But this is from 2023!!
water_table <- data.frame(read_excel(paste(project_root,"/metadata/environmental_data/Water_table_measurements.xlsx", sep = ""), sheet = 1), stringsAsFactors = F, check.names = F)

# remove empty rows
locs_na <- which(is.na(water_table$Päivämäärä))
water_table <- water_table[-locs_na,]

# add month
water_table$Month <- as.numeric(unlist(lapply(X = strsplit(x = as.character(water_table$Päivämäärä), split = "-"), function(x) x[2])))

# add day
water_table$Day <- as.numeric(unlist(lapply(X = strsplit(x = as.character(water_table$Päivämäärä), split = "-"), function(x) x[3])))

# add year
water_table$Year <- as.numeric(unlist(lapply(X = strsplit(x = as.character(water_table$Päivämäärä), split = "-"), function(x) x[1])))

# keep only year 2021, the year of sampling
water_table_2021 <- water_table[water_table$Year==2021,]

# rename finnish plot variable to english
colnames(water_table_2021)[2] <- "Plot"

# read-in prepared metadata
load(paste(project_root,"/metadata/Study_Metadata.RData", sep = ""))
rownames(metadata) <- metadata$Plot

# join together
water_table_joined <- merge(
  water_table_2021,
  metadata,
  by.x = "Plot",
  by.y = "Plot",
  all.x = TRUE
)

table(water_table_joined$Month)
#7   8   9  10 
#144 180 180 108 
# we have data starting from july

# separate into months, october not needed here
water_table_july <- water_table_joined[water_table_joined$Month==7,]
water_table_august <- water_table_joined[water_table_joined$Month==8,]
water_table_september <- water_table_joined[water_table_joined$Month==9,]

# plot 

# make a lineplot for each plot / site

# all the unique plots
uniq_plots <- unique(water_table_july$Plot)

# how many timepoints do we have altogether
nr_timepoints <- (length(which(water_table_july$Plot==uniq_plots[1])) +
                    length(which(water_table_august$Plot==uniq_plots[1])) +
                    length(which(water_table_september$Plot==uniq_plots[1])))

# # for each plot, gather a timeseries
plot_water_levels <- data.frame(matrix(nrow = length(uniq_plots), ncol = nr_timepoints))

for(i in 1:length(uniq_plots)){
  
  # get the july measurements
  a <- water_table_july[water_table_july$Plot==uniq_plots[i],]
  # order based on day
  a <- a[order(a$Day),]
  
  # get the august measurements
  b <- water_table_august[water_table_august$Plot==uniq_plots[i],]
  # order based on day
  b <- b[order(b$Day),]
  
  # get the september measurements
  c <- water_table_september[water_table_september$Plot==uniq_plots[i],]
  # order based on day
  c <- c[order(c$Day),]
  
  # combine into one variable
  plot_water_level <- c(as.numeric(a$`Water level from ground level`),
                        as.numeric(b$`Water level from ground level`),
                        as.numeric(c$`Water level from ground level`))
  
  # save into the common data frame for all plots
  plot_water_levels[i,] <- plot_water_level
}
# add rownames
rownames(plot_water_levels) <- uniq_plots

# add measurement dates as column names
colnames(plot_water_levels) <- c(paste(a$Day, a$Month, sep = "_"),
                                 paste(b$Day, b$Month, sep = "_"),
                                 paste(c$Day, c$Month, sep = "_"))

# separate by the exlusion treatment (here grazed and ungrazed)
all(uniq_plots==rownames(metadata)) # TRUE
grazed_plots <- which(metadata$Grazing=="grazed")
unrazed_plots <- which(metadata$Grazing=="ungrazed")

plot_water_levels_grazed <- plot_water_levels[grazed_plots,]
plot_water_levels_ungrazed <- plot_water_levels[unrazed_plots,]

# change into plotting directory
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

# plot line plots
{
  pdf("water_table_site_lineplots.pdf", width = 11.7, height = 8.3)
  
  # grazed
  plot(as.numeric(plot_water_levels_grazed[1,]), type="l",lwd=2, xlab = "Date", xaxt="n", main="", ylab = "Water level from ground level", cex.axis = 1.4, cex.lab = 1.5, ylim = c(-4.75,8), col="tan1")
  for(i in seq(from=2, to=nrow(plot_water_levels_grazed))){
    lines(as.numeric(plot_water_levels_grazed[i,]), lwd=2, col="tan1")
  }
  axis(side = 1, at = seq(1:14), labels = gsub("_", ".", colnames(plot_water_levels_grazed)), cex = 1.5, cex.axis = 1.4, cex.lab = 1.5)
  
  # ungrazed, add to the same plot
  for(i in seq(from=1, to=nrow(plot_water_levels_ungrazed))){
    lines(as.numeric(plot_water_levels_ungrazed[i,]), lwd=2, col="forestgreen")
  }
  dev.off()
}

# also plot mean trajectories

# calculate means, standard deviations and standard errors for each date
grazed   <- col_stats(plot_water_levels_grazed)
ungrazed <- col_stats(plot_water_levels_ungrazed)

# get stats for all months
months <- colnames(plot_water_levels_grazed)
x <- seq_along(months)

# calculte a range for the plot, use standard deviation to show variability for each date.
ylim <- range(
  grazed$mean - grazed$sd,
  grazed$mean + grazed$sd,
  ungrazed$mean - ungrazed$sd,
  ungrazed$mean + ungrazed$sd,
  na.rm = TRUE
)

# plot mean trajectories
{
  pdf("water_table_mean_trajectories.pdf", width = 11.7, height = 8.3)
  plot(
    x, grazed$mean,
    type = "n",
    ylim = ylim,
    xaxt = "n",
    xlab = "",
    cex.axis = 1.4, cex.lab = 1.5,
    ylab = "Water level from ground level",
    main = "Mean water level (± SD)"
  )
  
  # grazed SD ribbon + line
  polygon(
    c(x, rev(x)),
    c(grazed$mean - grazed$sd, rev(grazed$mean + grazed$sd)),
    col = adjustcolor("tan1", alpha.f = 0.3),
    border = NA
  )
  lines(x, grazed$mean, lwd = 3, col = "tan1")
  
  # Ungrazed SD ribbon + line
  polygon(
    c(x, rev(x)),
    c(ungrazed$mean - ungrazed$sd, rev(ungrazed$mean + ungrazed$sd)),
    col = adjustcolor("forestgreen", alpha.f = 0.3),
    border = NA
  )
  lines(x, ungrazed$mean, lwd = 3, col = "forestgreen")
  
  axis(1, at = x, labels = gsub("_", ".", months), las = 1, cex = 1.5, cex.axis = 1.4, cex.lab = 1.5)
  
  legend(
    "topleft",
    legend = c("Outside", "Inside"),
    col = c("tan1", "forestgreen"),
    lwd = 3,
    bty = "n",
    cex = 1.4
  )
  dev.off()
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures

