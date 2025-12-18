# a small script to plot Methanobrevibacter and Methanobacteriaceae in the metatranscriptomics data

# load libraries
library(phyloseq)
library(microbiome)
library(grid)
library(lmerTest)

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load metagenomics data
load(paste(project_root, "/metagenomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = ""))

# explore Methanobrevibacter
grep("Methanobrevibacter", tax_table_all$Genus)
# not detected

# load metatranscriptomics data
load(paste(project_root, "/metatranscriptomics/phyloflash/downstream/Otu_Tax_Tables_Parsed_Filtered.RData", sep = ""))

# explore Methanobrevibacter
locs_int <- grep("Methanobrevibacter", tax_table_all$Genus)
# it is detected in metatranscriptomics data

# make a phyloseq object - use the option where eukaryotes are summarized
taxa <- tax_table(as.matrix(tax_table_all))
otus <- otu_table(as.matrix(otu_table_all), taxa_are_rows = TRUE)
samples <- sample_data(metadata)
psq_all <- phyloseq(taxa, otus, samples)

# summarize to genus level
psq_genus <- aggregate_taxa(x = psq_all, level = "Genus")

# make compositional
psq_genus_rel <- microbiome::transform(x = psq_genus, transform = "compositional")

# clr transform
psq_genus_clr <- microbiome::transform(x = psq_genus_rel, transform = "clr")

# family level also
psq_family <- aggregate_taxa(x = psq_all, level = "Family")

# make compositional
psq_family_rel <- microbiome::transform(x = psq_family, transform = "compositional")

# clr transform
psq_family_clr <- microbiome::transform(x = psq_family_rel, transform = "clr")

# get Methanobrevibacter
tax_table_genus <- data.frame(tax_table(psq_genus_rel), stringsAsFactors = F, check.names = F)
otu_table_genus <- data.frame(otu_table(psq_genus_rel), stringsAsFactors = F, check.names = F)
otu_table_genus <- otu_table_genus * 100

# get also clr data
otu_table_genus_clr <- data.frame(otu_table(psq_genus_clr), stringsAsFactors = F, check.names = F)
locs_int <- grep("Methanobrevibacter", tax_table_genus$Genus)

# Methanobrevibacter summarized to genus level
methanobrevibacter_rel_sum <- otu_table_genus[locs_int,]
methanobrevibacter_clr_sum <- otu_table_genus_clr[locs_int,]

# Methanobacteriaceae
tax_table_family <- data.frame(tax_table(psq_family_rel), stringsAsFactors = F, check.names = F)
otu_table_family <- data.frame(otu_table(psq_family_rel), stringsAsFactors = F, check.names = F)
otu_table_family <- otu_table_family * 100

# get also clr data
otu_table_family_clr <- data.frame(otu_table(psq_family_clr), stringsAsFactors = F, check.names = F)
locs_int <- grep("Methanobacteriaceae", tax_table_family$Family)

# Methanobacteriaceae summarized to family level
methanobacteriaceae_rel_sum <- otu_table_family[locs_int,]
methanobacteriaceae_clr_sum <- otu_table_family_clr[locs_int,]

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# boxplot - plot one extra empty plot for easier postprocessing in inkscape

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))


# all relevant together for grazing
{
  pdf(file = "rumen_methanogens_boxplot_lmm.pdf", width = 11.7, height = 8.3, onefile = T)
  par(mar=c(10.1, 4.1, 4.1, 2.1))
  par(mfrow=c(2,5))
  plot(1,1)
  vals <- as.numeric(methanobrevibacter_rel_sum)
  g <- as.factor(meta_use$grazing)
  r_boxplot <- range(vals, na.rm = T)
  main_text <- "Methanobrevibacter"
  
  # italize
  main_text <- bquote(italic(.(main_text)))
  
  boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, ylab = "Relative abundance (%)", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  
  # use similar lmm models as elsewhere to get significance values
  vals <- as.numeric(methanobrevibacter_clr_sum)
  temp_data <- cbind(vals, meta_use)
  colnames(temp_data)[1] <- "vals"
  set.seed(1)
  mod_methanobrevibacter <- lmer(vals ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
  sum_temp <- summary(mod_methanobrevibacter)
  mtext(paste("P-value:", round(sum_temp$coefficients[2,5],3)), cex=0.8)
  
  vals <- as.numeric(methanobacteriaceae_rel_sum)
  g <- as.factor(meta_use$grazing)
  r_boxplot <- range(vals, na.rm = T)
  main_text <- "Methanobacteriaceae"
  
  # italize
  main_text <- bquote(italic(.(main_text)))
  
  boxplot(vals~g, col=c("forestgreen", "tan1"), xlab = "", ylim = r_boxplot, ylab = "Relative abundance (%)", main = main_text, cex.main = 0.7, cex.axis = 1.4, cex.lab = 1.5, outline = F, frame = FALSE)
  for(j in 1:length(unique(g))){
    points(x = jitter(rep(as.numeric(unique(g)[j]),length(which(g%in%unique(g)[j]))), factor = 2), y=vals[which(g%in%unique(g)[j])], pch=19, cex=0.7)
  }
  
  # use similar lmm models as elsewhere to get significance values
  vals <- as.numeric(methanobacteriaceae_clr_sum)
  temp_data <- cbind(vals, meta_use)
  colnames(temp_data)[1] <- "vals"
  set.seed(1)
  mod_methanobacteriaceae <- lmer(vals ~ grazing + treatment + (1|veg_clusters), REML=FALSE, data = temp_data)
  sum_temp <- summary(mod_methanobacteriaceae)
  mtext(paste("P-value:", round(sum_temp$coefficients[2,5],3)), cex=0.8)
  
  dev.off()
}

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures