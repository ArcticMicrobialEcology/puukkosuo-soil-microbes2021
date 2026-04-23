# a small script to plot the data from the soil moisture loggers from Puukkosuo from Summer 2021

# define some custom functions
# a function to get mean and sd
mean_sd <- function(x) {
  x <- as.numeric(x)
  c(
    mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE)
  )
}

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

moisture_data <- read.csv(paste(project_root, "/metadata/environmental_data/Oulanka_fen_soilmoisture_2021_rawdata_cleaned.csv", sep = ""))

# add month
moisture_data$month <- as.numeric(unlist(lapply(strsplit(x = moisture_data$TIMESTAMP, split = "-"), function(x) x[2])))

# add year
moisture_data$year <- as.numeric(unlist(lapply(strsplit(x = moisture_data$TIMESTAMP, split = "-"), function(x) x[1])))

# add date
# do in two parts
date <- unlist(lapply(strsplit(x = moisture_data$TIMESTAMP, split = "-"), function(x) x[3]))
date <- as.numeric(unlist(lapply(strsplit(x = date, split = "T"), function(x) x[1])))
moisture_data$date <- date

# separate to grazed and ungrazed
grazed_moisture <- moisture_data[which(moisture_data$grazing=="grazed"),]
table(grazed_moisture$plot_number)
table(grazed_moisture$month)
table(grazed_moisture$year)

ungrazed_moisture <- moisture_data[which(moisture_data$grazing=="ungrazed"),]
table(ungrazed_moisture$plot_number)
table(ungrazed_moisture$month)
table(ungrazed_moisture$year)

# only from 14 plots from the ungrazed plots or inside the fence, 4 plots had sensor malfunctions in the measurements

# plot

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))

{
  pdf(file = "soil_moisture_mean_trajectories.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(8.1, 4.1, 4.1, 2.1))
  # get mean trajectories for months
  months <- unique(moisture_data$month)
  
  # grazed
  grazed_stats <-data.frame(t(sapply(months, function(m) {
    mean_sd(grazed_moisture[grazed_moisture$month == m, 3])
  })), stringsAsFactors = FALSE, check.names = FALSE)
  
  # ungrazed
  ungrazed_stats <-data.frame(t(sapply(months, function(m) {
    mean_sd(ungrazed_moisture[ungrazed_moisture$month == m, 3])
  })), stringsAsFactors = FALSE, check.names = FALSE)
  
  # calculate p-value with wilcoxon text, this now compares all the moisture data values for a month between outside and inside the exclosure
  p_vals <- sapply(months, function(m){
    wilcox.test(x = as.numeric(grazed_moisture[grazed_moisture$month == m, 3]), 
                y=as.numeric(ungrazed_moisture[ungrazed_moisture$month == m, 3]), 
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
  
  plot(
    x, grazed_stats[, "mean"],
    type = "l", lwd = 3, col = "tan1",
    ylim = y_lim,
    xaxt = "n",
    xlab = "", ylab = "Soil moisture mV",
    main = "",
    cex.axis = 1.4, cex.lab = 1.5
  )
  
  # SE ribbon (grazed)
  polygon(
    c(x, rev(x)),
    c(grazed_stats[, "mean"] - grazed_stats[, "sd"],
      rev(grazed_stats[, "mean"] + grazed_stats[, "sd"])),
    col = adjustcolor("tan1", alpha.f = 0.3),
    border = NA
  )
  
  # Ungrazed line
  lines(x, ungrazed_stats[, "mean"], lwd = 3, col = "forestgreen")
  
  # SE ribbon (ungrazed)
  polygon(
    c(x, rev(x)),
    c(ungrazed_stats[, "mean"] - ungrazed_stats[, "sd"],
      rev(ungrazed_stats[, "mean"] + ungrazed_stats[, "sd"])),
    col = adjustcolor("forestgreen", alpha.f = 0.3),
    border = NA
  )
  
  # add months axis
  axis(1, at = x, labels = c("May", "June", "July", "August", "September"), las = 2, cex.axis = 1.4, cex.lab = 1.5)
  legend(
    "topright",
    legend = c("Outside", "Inside"),
    col = c("tan1", "forestgreen"),
    lwd = 3, bty = "n",
    cex = 1.5
  )
  
  # add p-vals as axis above the figure
  p_vals <- as.character(round(p_vals, 3))
  if(any(p_vals=="0")){p_vals[which(p_vals=="0")] <- "<0.001"}
  axis(3, at = x, labels = p_vals, las = 2, cex.axis = 1.4, cex.lab = 1.5)
  
  dev.off()
}

# plot lines for each plot for the whole summer aswell

# grazed
unique_plots <- unique(grazed_moisture$plot_number)
grazed_df <- matrix(nrow = length(unique_plots), ncol = as.numeric(table(grazed_moisture$plot_number)[1])) # same amount of measurements for each plot

# convert the timestamp to POSIXct
grazed_moisture$TIMESTAMP_pos <- as.POSIXct(grazed_moisture$TIMESTAMP, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

# gather the data
for(i in 1:length(unique_plots)){
  # pick all values for the plot
  temp <- grazed_moisture[which(grazed_moisture$plot_number==unique_plots[i]),]
  
  if(i==1){
    # order based on timestamp
    temp <- temp[order(temp$TIMESTAMP_pos, decreasing = FALSE), ]
    colnames(grazed_df) <- temp$TIMESTAMP
  }else{
    # put in the same order
    rownames(temp) <- temp$TIMESTAMP
    temp <- temp[colnames(grazed_df),]
  }
  
  # save into the common df
  grazed_df[i,] <- as.numeric(temp$soilmoist)
}
rownames(grazed_df) <- unique_plots
#colnames(grazed_df) <- paste(temp$date, temp$month, sep = ".")

# ungrazed
unique_plots <- unique(ungrazed_moisture$plot_number)
ungrazed_df <- matrix(nrow = length(unique_plots), ncol = as.numeric(table(ungrazed_moisture$plot_number)[1])) # same amount of measurements for each plot

# convert the timestamp to POSIXct
ungrazed_moisture$TIMESTAMP_pos <- as.POSIXct(ungrazed_moisture$TIMESTAMP, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

for(i in 1:length(unique_plots)){
  # pick all values for the plot
  temp <- ungrazed_moisture[which(ungrazed_moisture$plot_number==unique_plots[i]),]
  
  if(i==1){
    # order based on timestamp
    temp <- temp[order(temp$TIMESTAMP_pos, decreasing = FALSE), ]
    colnames(ungrazed_df) <- temp$TIMESTAMP
  }else{
    # put in the same order
    rownames(temp) <- temp$TIMESTAMP
    temp <- temp[colnames(ungrazed_df),]
  }
  
  # save into the common df
  ungrazed_df[i,] <- as.numeric(temp$soilmoist)
}
rownames(ungrazed_df) <- unique_plots
# colnames(ungrazed_df) <- paste(temp$date, temp$month, sep = ".")

# check that timestamps match
all(colnames(grazed_df)==colnames(ungrazed_df)) # TRUE
# moisture data for all samples in the same order

# plot
{
  pdf(file = "soil_moisture_individual_plot_trajectories.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(8.1, 4.1, 4.1, 2.1))
  
  y_lim <- range(c(unlist(grazed_df), unlist(ungrazed_df)), na.rm = T)
  plot(as.numeric(grazed_df[1,]), type="l", lwd=2, xlab = "", xaxt="n", main="", ylab = "Soil moisture mV", ylim = y_lim, col="tan1", cex.axis = 1.4, cex.lab = 1.5)
  for(i in seq(from=2, to=nrow(grazed_df))){
    lines(as.numeric(grazed_df[i,]), lwd=2, col="tan1")
  }
  
  # add month as x axis for the first measurements for each month
  axis(side = 1, at = c(0,which(diff(temp$month)==1)), labels = c("May", "June", "July", "August", "September"), las=2, cex.axis = 1.4, cex.lab = 1.5)
  
  
  # ungrazed
  for(i in seq(from=1, to=nrow(ungrazed_df))){
    lines(as.numeric(ungrazed_df[i,]), lwd=2, col="forestgreen")
  }
  
  legend(
    "topright",
    legend = c("Outside", "Inside"),
    col = c("tan1", "forestgreen"),
    lwd = 3,
    bty = "n",
    cex = 1.5
  )
  dev.off()

}

# plot mean trajectories still for each timestamp
# calculate means, standard deviations and standard errors for each date
grazed   <- col_stats(grazed_df)
ungrazed <- col_stats(ungrazed_df)

# get stats for all months
timepoints <- colnames(grazed_df)
x <- seq_along(timepoints)

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
  pdf("soil_moisture_mean_trajectories_over_all_timepoints.pdf", width = 11.7, height = 8.3)
  par(mar=c(8.1, 4.1, 4.1, 2.1))
  plot(
    x, grazed$mean,
    type = "n",
    ylim = ylim,
    xaxt = "n",
    xlab = "",
    cex.axis = 1.4, cex.lab = 1.5,
    ylab = "Soil moisture mV",
    main = ""
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
  
  # add month as x axis for the first measurements for each month
  axis(side = 1, at = c(0,which(diff(temp$month)==1)), labels = c("May", "June", "July", "August", "September"), las=2, cex.axis = 1.4, cex.lab = 1.5)
  
  # axis(1, at = x, labels = gsub("_", ".", months), las = 1, cex = 1.5, cex.axis = 1.4, cex.lab = 1.5)
  
  legend(
    "topright",
    legend = c("Outside", "Inside"),
    col = c("tan1", "forestgreen"),
    lwd = 3,
    bty = "n",
    cex = 1.5
  )
  dev.off()
}
