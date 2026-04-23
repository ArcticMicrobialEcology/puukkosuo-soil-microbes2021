# a small script to produce some summary statistics for the MAGs

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load data for the MAGs
load("/scratch/project_2009164/2_OULANKA/Tommi/final/metagenomics/mag_based/final_mags/Final_MAG_Taxonomy_Abundance.RData")
load(paste(project_root, "/metagenomics/mag_based/final_mags/Final_MAG_Taxonomy_Mark_Gene.RData", sep = ""))

# produce a table and some summaries of MAG completeness etc for the manuscript.

# read-in thje sumamries of the bins from the different conditions
# grazed
setwd("/scratch/project_2009164/2_OULANKA/Tommi/final/metagenomics/mag_based/Grazed/MAGs/metabat2/summary_refined_mags")
grazed_bins <- read.csv(paste(project_root, "/metagenomics/mag_based/Grazed/MAGs/metabat2/summary_refined_mags/bins_summary.txt", sep = ""), sep = "\t")
rownames(grazed_bins) <- grazed_bins$bins

# match the order of the final GTDB annotated report
grazed_mags <- dereplication_report[which(dereplication_report$condition=="Grazed"),]
rownames(grazed_mags) <- grazed_mags$bin_name_in_condition
grazed_bins <- grazed_bins[rownames(grazed_mags),]

# ungrazed
ungrazed_bins <- read.csv(paste(project_root, "/metagenomics/mag_based/UnGrazed/MAGs/metabat2/summary_refined_mags/bins_summary.txt", sep = ""), sep = "\t")
rownames(ungrazed_bins) <- ungrazed_bins$bins

# match the order of the final GTDB annotated report
ungrazed_mags <- dereplication_report[which(dereplication_report$condition=="UnGrazed"),]
rownames(ungrazed_mags) <- ungrazed_mags$bin_name_in_condition
ungrazed_bins <- ungrazed_bins[rownames(ungrazed_mags),]

# combine and parse
all_bins <- data.frame(matrix(nrow = nrow(dereplication_report), ncol = ncol(grazed_bins))) 
colnames(all_bins) <- colnames(grazed_bins)

all(dereplication_report$bin_name_in_condition[which(dereplication_report$condition=="Grazed")]==rownames(grazed_bins)) #TRUE
all_bins[which(dereplication_report$condition=="Grazed"), ] <- grazed_bins

all(dereplication_report$bin_name_in_condition[which(dereplication_report$condition=="UnGrazed")]==rownames(ungrazed_bins)) #TRUE
all_bins[which(dereplication_report$condition=="UnGrazed"), ] <- ungrazed_bins
rownames(all_bins) <- rownames(dereplication_report)
all_bins <- cbind(all_bins, dereplication_report)

# load and add taxonomy and MAG information as rownames
load("/scratch/project_2009164/2_OULANKA/Tommi/final/manuscript/fig6_related/modified_gtdb_tax_silva.RData")
load(paste(project_root, "/downstream/manuscript_figures/figure6/modified_gtdb_tax_silva.RData", sep = ""))
all(all_bins$cluster==rownames(gtdbtk_report)) # TRUE

m_num <- seq(1:nrow(gtdbtk_report))
r_names <- paste("MAG",m_num, "; ",gtdbtk_report$Phylum, "; ", gtdbtk_report$Order, "; ", gtdbtk_report$Genus, sep = "")
rownames(all_bins) <- r_names

# save the table

# change directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/supplementary", sep = ""))
save(all_bins, file = "All_MAGs_info.RData")
write.csv(x = all_bins, file = "All_MAGs_info.csv", row.names = T, col.names = T, quote = F)

# make also a summary table
mag_summary_table <- rbind(summary(all_bins$total_length),
                           summary(all_bins$num_contigs),
                           summary(all_bins$N50),
                           summary(all_bins$GC_content),
                           summary(all_bins$percent_completion),
                           summary(all_bins$percent_redundancy))
rownames(mag_summary_table) <- c("total_length", "num_contigs", "N50", "GC_content", "percent_completion", "percent_redundancy")                           
colnames(mag_summary_table) <- names(summary(all_bins$total_length))
save(mag_summary_table, file = "MAG_summary.RData")
write.csv(x = mag_summary_table, file = "MAG_summary.csv", row.names = T, col.names = T, quote = F)

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess as needed