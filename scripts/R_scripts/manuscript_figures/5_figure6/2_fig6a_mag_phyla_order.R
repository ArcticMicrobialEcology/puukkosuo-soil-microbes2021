# a small script to plot the most abundant taxonomic annotations of  the MAGs

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(phyloseq)
library(microbiome)
library(grid)
library(RColorBrewer)

# load the processed taxonomy data for the MAGs
load(paste(project_root, "/metagenomics/mag_based/final_mags/Final_MAG_Taxonomy_Mark_Gene.RData", sep = ""))

# phylum level
# for plotting, collapse Desulfobacterota_B etc. into Desulfobacterota (A-E) and Myxococcota_ into Myxococcota 
# to keep annotation similar to the read-based analysis
all_phylum <- gtdbtk_report$Phylum
all_phylum[grep("Desulfobacterota", all_phylum)] <- "Desulfobacterota (A-E)"
all_phylum[grep("Myxococcota", all_phylum)] <- "Myxococcota"
taxa_phylum <- sort(table(all_phylum), decreasing = T)
length(taxa_phylum) #13

# load plot colors - use the same coloring scheme as in figure 3
load(paste(project_root, "/downstream/manuscript_figures/figure3/phylum_plot_colors.RData", sep = ""))
# load more colors
load(paste(project_root, "/downstream/manuscript_figures/figure6/new_colors.RData", sep = ""))

# plot phylum first
# set colors
plot_cols <- new20[1:length(taxa_phylum)]
names(plot_cols) <- names(taxa_phylum)
plot_cols[which(names(plot_cols)%in%names(cols_use))]<- cols_use[names(plot_cols)[which(names(plot_cols)%in%names(cols_use))]]

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure6", sep = ""))

# plot phylum level
{
  pdf(file = "barplot_MAG_taxa_phylum_MG.pdf", width = 11.7, height = 8.3)
  
  # desired bar width in inches
  desired_bar_width_inches <- 1
  
  # total page width in inches
  total_page_width_inches <- 11.7
  
  # calculate the left and right margins to center the bar
  margin_inches <- (total_page_width_inches - desired_bar_width_inches) / 2
  
  # set the plot dimensions in inches
  options(repr.plot.width = total_page_width_inches, repr.plot.height = 4, repr.plot.units = "in") 
  
  # calculate relative bar width and margins
  bar_width_relative <- desired_bar_width_inches / total_page_width_inches
  margin_relative <- margin_inches / total_page_width_inches
  
  # plot
  barplot(as.matrix(taxa_phylum), main = "" , width = bar_width_relative, col=plot_cols,space = c(margin_relative, 0), xlim = c(0,1), border="white", las=2, xlab="", ylab = "Number of MAGs", ylim=c(0,120))
  
  # reset graphics options
  options(repr.plot.width = NULL, repr.plot.height = NULL, repr.plot.units = NULL)
  
  # legends for next page
  grid.newpage()
  legend("top", legend = lapply(names(plot_cols), function(x) bquote(italic(.(x)))), fill = plot_cols, ncol=3)
  dev.off()
}

# plot order level
# add parsed phyla to taxonomy table
gtdbtk_report$Parsed_Phylum <- all_phylum
save(gtdbtk_report, file = "modified_gtdb_tax_silva.RData")

# get order level taxa
taxa_order <- sort(table(gtdbtk_report$Order), decreasing = T)
length(taxa_order) #37

# truncate to those that appear more than once
taxa_order <- taxa_order[which(taxa_order>1)]
length(taxa_order) # 20
sum(taxa_order) #96

# add other category - number of MAGs belonging to this category is 17
taxa_order <- c(taxa_order,17)
names(taxa_order)[length(taxa_order)] <- "Other"

# load plot colors - use the same coloring scheme as in figure 3
load(paste(project_root, "/downstream/manuscript_figures/figure3/order_plot_colors.RData", sep = ""))

# add one more color
new20 <- c(new20, "gray8")

# plot order
plot_cols <- new20[1:length(taxa_order)]
names(plot_cols) <- names(taxa_order)
plot_cols[which(names(plot_cols)%in%names(cols_use))]<- cols_use[names(plot_cols)[which(names(plot_cols)%in%names(cols_use))]]

# add phylum names
phyl_names <- character(length(plot_cols))
for(i in 1:(length(plot_cols)-1)){
  locs_int <- which(gtdbtk_report$Order%in%names(plot_cols)[i])
  phyl_names[i] <- unique(gtdbtk_report$Parsed_Phylum[locs_int])
}
phyl_names[length(phyl_names)] <- "Other"
names(plot_cols) <- paste(phyl_names, names(plot_cols), sep = "; ")

# plot
{
  pdf(file = "barplot_MAG_taxa_order_MG.pdf", width = 11.7, height = 8.3)
  
  # desired bar width in inches
  desired_bar_width_inches <- 1
  
  # total page width in inches
  total_page_width_inches <- 11.7
  
  # calculate the left and right margins to center the bar
  margin_inches <- (total_page_width_inches - desired_bar_width_inches) / 2
  
  # set the plot dimensions in inches
  options(repr.plot.width = total_page_width_inches, repr.plot.height = 4, repr.plot.units = "in") 
  
  # calculate relative bar width and margins
  bar_width_relative <- desired_bar_width_inches / total_page_width_inches
  margin_relative <- margin_inches / total_page_width_inches
  
  # plot
  barplot(as.matrix(taxa_order), main = "" , width = bar_width_relative, col=plot_cols,space = c(margin_relative, 0), xlim = c(0,1), border="white", las=2, xlab="", ylab = "Number of MAGs", ylim=c(0,120))
  
  # reset graphics options
  options(repr.plot.width = NULL, repr.plot.height = NULL, repr.plot.units = NULL)
  
  # legends for next page
  grid.newpage()
  legend("top", legend = lapply(names(plot_cols), function(x) bquote(italic(.(x)))), fill = plot_cols, ncol=3)
  dev.off()
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures
