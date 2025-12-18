# a small script to plot the most abundant and expressed taxa in the MG and MT datas

# load libraries
library(phyloseq)
library(microbiome)
library(grid)

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load datas
# metagenomics
load(paste(project_root, "/metagenomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = ""))

# here we have compositional data already aggregated to the order level
psq_order_mg <- psq_order

# compose order level data with proper taxonomy, as aggregate_rare by microbiome does not preserve taxonomy tables well
psq_order_mg_tax <- aggregate_taxa(x = psq, level = "Order")

# get phylum levels as well
psq_phylum_mg <- psq_phylum

# metatranscriptomics 
load(paste(project_root, "/metatranscriptomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = ""))

# here we have compositional data already aggregated to the order level
psq_order_mt <- psq_order

# compose order level data with proper taxonomy, as aggregate_rare by microbiome does not preserve taxonomy tables well
psq_order_mt_tax <- aggregate_taxa(x = psq, level = "Order")

# get phylum levels as well
psq_phylum_mt <- psq_phylum

### start plotting
## order level
# get compositional datas
# metagenomics
data_order_mg <- data.frame(otu_table(psq_order_mg), stringsAsFactors = F, check.names = F) * 100

# metatranscriptomics
data_order_mt <- data.frame(otu_table(psq_order_mt), stringsAsFactors = F, check.names = F) * 100

# get 10 most abundant taxa
# metagenomics
data_temp_mg <- data_order_mg
unkn_loc <- grep("unknown|Other", rownames(data_temp_mg))
if(length(unkn_loc)>0){data_temp_mg <- data_temp_mg[-unkn_loc,]}
top_taxa_mg <- names(sort(rowSums(data_temp_mg), decreasing = T))[1:10]

# metatranscriptomics
data_temp_mt <- data_order_mt
unkn_loc <- grep("unknown|Other", rownames(data_temp_mt))
if(length(unkn_loc)>0){data_temp_mt <- data_temp_mt[-unkn_loc,]}
top_taxa_mt <- names(sort(rowSums(data_temp_mt), decreasing = T))[1:10]

# get all top taxa over both MG and MT
all_top_taxa <- unique(c(top_taxa_mg, top_taxa_mt))

# define some divergent colors for plotting
plot_colors <- c("purple", "cornflowerblue","mistyrose2","darkgrey","deeppink","brown","violet","firebrick1","chocolate3","forestgreen","darkslateblue","darkolivegreen",
                 "chocolate1", "cadetblue1", "darkgoldenrod3", "blue", "darkolivegreen3", "yellow", "burlywood2", "black", "green")

cols_use <- plot_colors[1:length(all_top_taxa)]
names(cols_use) <- all_top_taxa

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure3", sep = ""))

# save the used colors for later use
save(plot_colors, cols_use, file = "order_plot_colors.RData")

# plot - using similar color scheme for metagenomics and metatranscriptomics
# metagenomics
data_plot <- data_order_mg[top_taxa_mg,]

# add phylum names to order data
tax_mg <- data.frame(tax_table(psq_order_mg_tax), stringsAsFactors = F, check.names = F)
tax_mg <- tax_mg[rownames(data_plot),]

# map Silva phylum annotations to GTDB phylum annotations for the visualization
# load the precompiled mapping file
load("Silva_GTDB_phylum_mapping.RData")
all_phyla <- tax_mg$Phylum
rownames(silva_gtdb_map) <- silva_gtdb_map$SILVA
all_phyla <- silva_gtdb_map[all_phyla, 2]
namn <- paste(rownames(data_plot), all_phyla, sep = "; ")

# plot
{
  pdf(file = "barplot_top10_taxa_order_MG.pdf", width = 11.7, height = 8.3)
  par(mar=c(6.1, 4.1, 4.1, 2.1))
  barplot(as.matrix(data_plot), col=cols_use[rownames(data_plot)], ylim=c(0,40), ylab = "Abundance (%)", main = "", las=2)
  grid.newpage()
  legend("top", legend = lapply(namn, function(x) bquote(italic(.(x)))), fill = cols_use[rownames(data_plot)], ncol=3)
  dev.off()
}

# metatranscriptomics
data_plot <- data_order_mt[top_taxa_mt,]

# add phylum names to order data
tax_mt <- data.frame(tax_table(psq_order_mt_tax), stringsAsFactors = F, check.names = F)
tax_mt <- tax_mt[rownames(data_plot),]

# map Silva phylum annotations to GTDB phylum annotations for the visualization
all_phyla <- tax_mt$Phylum
all_phyla <- silva_gtdb_map[all_phyla, 2]
namn <- paste(rownames(data_plot), all_phyla, sep = "; ")

# plot
{
  pdf(file = "barplot_top10_taxa_order_MT.pdf", width = 11.7, height = 8.3)
  par(mar=c(6.1, 4.1, 4.1, 2.1))
  barplot(as.matrix(data_plot), col=cols_use[rownames(data_plot)], ylim=c(0,40), ylab = "Abundance (%)", main = "", las=2)
  grid.newpage()
  legend("top", legend = lapply(namn, function(x) bquote(italic(.(x)))), fill = cols_use[rownames(data_plot)], ncol=3)
  dev.off()
}

## phylum level
# get compositional data
data_phylum_mg <- data.frame(otu_table(psq_phylum_mg), stringsAsFactors = F, check.names = F) * 100
data_phylum_mt <- data.frame(otu_table(psq_phylum_mt), stringsAsFactors = F, check.names = F) * 100

# get 10 most abundant taxa in samples
# metagenomics
data_temp_mg <- data_phylum_mg
unkn_loc <- grep("unknown|Other", rownames(data_temp_mg))
if(length(unkn_loc)>0){data_temp_mg <- data_temp_mg[-unkn_loc,]}
top_taxa_mg <- names(sort(rowSums(data_temp_mg), decreasing = T))[1:10]

# metatranscriptomics
data_temp_mt <- data_phylum_mt
unkn_loc <- grep("unknown|Other", rownames(data_temp_mt))
if(length(unkn_loc)>0){data_temp_mt <- data_temp_mt[-unkn_loc,]}
top_taxa_mt <- names(sort(rowSums(data_temp_mt), decreasing = T))[1:10]

# get all top taxa over both MG and MT
all_top_taxa <- unique(c(top_taxa_mg, top_taxa_mt))

# define the colors to be used
cols_use <- plot_colors[1:length(all_top_taxa)]
names(cols_use) <- all_top_taxa

# plot metagenomics
data_plot <- data_phylum_mg[top_taxa_mg,]

# map to GTDB phyla
namn <-silva_gtdb_map[rownames(data_plot),2]
{
  pdf(file = "barplot_top10_taxa_phylum_MG.pdf", width = 11.7, height = 8.3)
  par(mar=c(6.1, 4.1, 4.1, 2.1))
  barplot(as.matrix(data_plot), col=cols_use[rownames(data_plot)], ylab = "Abundance (%)", main = "", las=2)
  grid.newpage()
  legend("top", legend = lapply(namn, function(x) bquote(italic(.(x)))), fill = cols_use[rownames(data_plot)], ncol=3)
  dev.off()
}

# metatranscriptomics
data_plot <- data_phylum_mt[top_taxa_mt,]

# map to GTDB phyla
namn <-silva_gtdb_map[rownames(data_plot),2]
{
  pdf(file = "barplot_top10_taxa_phylum_MT.pdf", width = 11.7, height = 8.3)
  par(mar=c(6.1, 4.1, 4.1, 2.1))
  barplot(as.matrix(data_plot), col=cols_use[rownames(data_plot)], ylab = "Abundance (%)", main = "", las=2)
  grid.newpage()
  legend("top", legend = lapply(namn, function(x) bquote(italic(.(x)))), fill = cols_use[rownames(data_plot)], ncol=3)
  dev.off()
}

# save the used colors but with GTDB phylum names - to be used also for MAG plotting
names(cols_use) <- silva_gtdb_map[names(cols_use),2]
save(plot_colors, cols_use, file = "phylum_plot_colors.RData")

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures