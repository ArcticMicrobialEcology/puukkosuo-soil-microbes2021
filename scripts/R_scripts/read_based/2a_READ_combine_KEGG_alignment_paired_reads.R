args <- commandArgs(trailingOnly = TRUE)
# requires the KEGG hits for R1 and the KEGG hits for R2 as the first and second input parameter 
# the third parameter needs to be the (downloaded) mapping file (tab delimited) between KEGG genes and KO groups (only two columns, gene and ko)
# as fourth parameter if the KEGG hits should be filtered (1) or not (0) with the defined default limits.
# currently the default limits: e-value 10^-6 and bitscore 50, change if needed
# fifth parameter needs to be the output directory (needs to exist) where the subdirectories for the results are created (if they don't exist) and results are saved

# e.g.
# args <- character(6)
# args[1] <- "A039-MB-M-TGACGGCCGT-AGTCAACCAT-Velmala-0121-A_kegg_R1.txt"
# args[2] <- "A039-MB-M-TGACGGCCGT-AGTCAACCAT-Velmala-0121-A_kegg_R2.txt"
# args[3] <- "/scratch/project_2007998/AGROBIO_Tommi/metadata/2025_01_21_KO_Mapping_Genes_Prokaryotes.txt"
# args[4] <- "0"
# args[5] <- "/scratch/project_2007998/AGROBIO_Tommi/metagenomics/kegg_diamond"

# a script to combine the results of the DIAMOND blastp KEGG homology search for the read pairs of paired reads

# load libraries
library(grr)
library(stringdist)

# read-in the input files, this format is fixed and follows the one defined in the shell script for the DIAMOND search run.

# DIAMOND results for R1
r1_func <- read.csv(as.character(args[1]), sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F, 
                    colClasses = c("character","character","numeric", "integer", "integer", "integer", 
                                   "integer", "integer","numeric","numeric", "numeric", "numeric", "integer"))
colnames(r1_func) <- c("qseqid", "stitle", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp", "slen")

# DIAMOND results for R2
r2_func <- read.csv(as.character(args[2]), sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F, 
                    colClasses = c("character","character","numeric", "integer", "integer", "integer", 
                                   "integer", "integer","numeric","numeric", "numeric", "numeric", "integer"))
colnames(r2_func) <- c("qseqid", "stitle", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp", "slen")

# what was imported?
print(paste("Read-in R1 file with", nrow(r1_func), "hits"))
print(paste("Read-in R2 file with", nrow(r2_func), "hits"))

# filter if desired, change the limits if needed
if(as.numeric(args[4]) == 1){
  r1_func <- r1_func[which(r1_func$evalue <= 10^(-6) & r1_func$bitscore >= 50),]
  r2_func <- r2_func[which(r2_func$evalue <= 10^(-6) & r2_func$bitscore >= 50),]
  
  print(paste("R1 hits remaining after filtering:", nrow(r1_func)))
  print(paste("R2 hits remaining after filtering:", nrow(r2_func)))
}

# get sample name - the naming of the files must be as defined in the shell script
sample_namn <- gsub("_kegg_R1.txt", "", args[1])

# parse gene level information.
print("Parsing input data")
r1_genes <- unlist(lapply(X = r1_func[,2], FUN =  function(x) strsplit(x = x[1], split = " ")[[1]][1]))
r2_genes <- unlist(lapply(X = r2_func[,2], FUN =  function(x) strsplit(x = x[1], split = " ")[[1]][1]))

# for each read pair, we want to calculate hits to exactly the same gene (or protein) only once. So if both read pairs have a hit in the same gene, remove pair2 hit.

# join gene information with other
r1_func <- cbind(r1_genes, r1_func)
r2_func <- cbind(r2_genes, r2_func)

colnames(r1_func)[1:2] <- c("gene", "read_location")
colnames(r2_func)[1:2] <- c("gene", "read_location")

# combine the hits from paired reads
print("Combining the hits from the paired reads")

# reads from the same location - paired reads
int_reads <- intersect(r1_func$read_location, r2_func$read_location)
rownames(r1_func) <- r1_func$read_location
rownames(r2_func) <- r2_func$read_location

# pick the paired reads
r1_func_int <- r1_func[int_reads,]
r2_func_int <- r2_func[int_reads,]

# R2 reads not having a hit for the R1 pair
r2_uniq <- r2_func[-which(r2_func$read_location%in%int_reads),]

# put reads from the same location for R1 and R2 in the same order
r1_func_int <- r1_func_int[order(r1_func_int$read_location),]
r2_func_int <- r2_func_int[order(r2_func_int$read_location),]

# this works, because there is only one (best) hit for each read and they are in the same order. So works with only one hit / read. 
locs_rem <- which(r2_func_int$gene==r1_func_int$gene)

# remove these duplicated gene hits as we don't want to count them
if(length(locs_rem)>0){
  r1_keep <- r1_func_int[-locs_rem,]
  r2_keep <- r2_func_int[-locs_rem,]
}else{
  r1_keep <- r1_func_int
  r2_keep <- r2_func_int
}

# now further filter the R2 reads to contain only those hits with KEGG or KO annotation differing from the KO annotation of the corresponding R1 pair

# read-in the mapping between KEGG genes and KO groups
PROKARYOTES.DAT <- read.csv(args[3], sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F)
colnames(PROKARYOTES.DAT)=c("gene", "KO")

# map the genes to ko groups in the filtered reads
# assign KOs, get results in a list as genes may have multiple KO annotations

# R1
kos1 <- matches(r1_keep$gene, PROKARYOTES.DAT$gene, all.x = TRUE, all.y = FALSE, list = TRUE, indexes = TRUE)
kos1 <- lapply(kos1, function(x) PROKARYOTES.DAT$KO[x])
names(kos1) <- r1_keep$read_location

# R2
kos2 <- matches(r2_keep$gene, PROKARYOTES.DAT$gene, all.x = TRUE, all.y = FALSE, list = TRUE, indexes = TRUE)
kos2 <- lapply(kos2, function(x) PROKARYOTES.DAT$KO[x])
names(kos2) <- r2_keep$read_location

# add also KO annotations to the larger dataframes for easier parsing later

# R1
all_kos1 <- matches(r1_func$gene, PROKARYOTES.DAT$gene, all.x = TRUE, all.y = FALSE, list = TRUE, indexes = TRUE)
all_kos1 <- lapply(all_kos1, function(x) PROKARYOTES.DAT$KO[x])
names(all_kos1) <- r1_func$read_location
all_kos1 <- unlist(lapply(all_kos1, function(x) paste(x, collapse = ";")))
r1_func$ko <- all_kos1

# R2, filtered hits
keep_kos2 <- unlist(lapply(kos2, function(x) paste(x, collapse = ";")))
r2_keep$ko <- keep_kos2

# R2, unique hits
unique_kos2 <- matches(r2_uniq$gene, PROKARYOTES.DAT$gene, all.x = TRUE, all.y = FALSE, list = TRUE, indexes = TRUE)
unique_kos2 <- lapply(unique_kos2, function(x) PROKARYOTES.DAT$KO[x])
names(unique_kos2) <- r2_uniq$read_location
unique_kos2 <- unlist(lapply(unique_kos2, function(x) paste(x, collapse = ";")))
r2_uniq$ko <- unique_kos2

# which R2 hits, don't have KO annotations, further examine these
l_r2 <- lengths(kos2)
expl_r2 <- which(l_r2==0) 
# there is most likely always those R2 hits with no KO annotations, no if sentence necessary?
# or in that case, this script will just break.

r1_expl <- r1_keep[expl_r2,]
r2_expl <- r2_keep[expl_r2,]

# these hits don't have KO annotation for R2 but their gene name for the paired reads is different than for corresponding R1 pair.
# Thus, they have the potential to be hits to truly different proteins (instead of hitting a very similar homolog than R1)
# following, explore the functionalities for the hits for these paired reads more closely.

# further filter these results based on delivered KEGG function similarity
r1_function <- unlist(lapply(X = r1_expl[,3], FUN =  function(x) {
  y <- strsplit(x = x[1], split = " ")[[1]]
  paste(y[3:length(y)], collapse = " ")
}))
r2_function <- unlist(lapply(X = r2_expl[,3], FUN =  function(x) {
  y <- strsplit(x = x[1], split = " ")[[1]]
  paste(y[3:length(y)], collapse = " ")
}))

# don't consider those hits in R2 with very similar function for the same readpair, as the paired reads are most likely just hitting very close homologs and these hits shouldn't be counted twice
func_similarity <- stringsim(a = r1_function, b = r2_function)
# osa similarity, the default method is used.
# The Levenshtein distance (method='lv') counts the number of deletions, insertions and substitutions necessary to turn b into a. 
# This method is equivalent to R's native adist function.
# The Optimal String Alignment distance (method='osa') is like the Levenshtein distance but also allows transposition of adjacent 
# characters. Here, each substring may be edited only once. (For example, a character cannot be transposed twice to move it forward in the string).

# remove those hits whos functionality string is more than 0.8 similar, ranges from 0 to 1
locs_rem <- which(func_similarity>0.8) 

if(length(locs_rem)>0){
  
  # the function of these are not known, keep these even if the function string are similar
  if(length(grep("hypothetical protein", r1_function[locs_rem], ignore.case = TRUE))>0){
    locs_rem <- locs_rem[-grep("hypothetical protein", r1_function[locs_rem], ignore.case = TRUE)] 
  }
  
  r1_expl <- r1_expl[-locs_rem,]
  r2_expl <- r2_expl[-locs_rem,]
}

# keep these hits for R2
all_keep <- c(rownames(r2_expl))

# remove these from the lists for further exploration
kos1_filt <- kos1[-expl_r2]
kos2_filt <- kos2[-expl_r2]

# all the remaining in the list to be investigated have KO annotation for R2
# which of these don't have annotation for R1? Keep the R2 annotations for all of these.
l_r1 <- lengths(kos1_filt)
keep_r2 <- which(l_r1==0)

# keep all of these
if(length(keep_r2)>0){
  all_keep <- c(all_keep, names(keep_r2))
  
  # remove these from the lists for further comparison
  kos1_filt <- kos1_filt[-keep_r2]
  kos2_filt <- kos2_filt[-keep_r2]
}
# now all the remaining hits have KO annotations for R1 and R2

# check which R2 KO annotations are truly different from the R1 annotations - keep only those for R2
# another option would be to keep all annotations but that would result in possibly different annotations for the same R1 genes - something we don't perhaps want.
found2 <- logical(length(kos2_filt))
for(i in 1:length(kos2_filt)){
  found2[i] <- any(kos2_filt[[i]]%in%kos1_filt[[i]])
}
names(found2) <- names(kos2_filt)
keep_r2 <- which(!found2)

all_keep <- c(all_keep, names(keep_r2)) 
# here we now have the names for all the R2 read hits we should keep

if(length(all_keep)>0){ 
  r2_keep <- r2_keep[all_keep,]
  colnames(r2_keep)[1] <- "gene"
  r2_func_filt <- rbind(r2_uniq, r2_keep)
} else { # we don't have any additional R2 hits we want to keep in addition to read pairs having hits only in R2
  r2_func_filt <- r2_uniq
}
# here we now have the R2 hits that either don't have a hit in R1 for the same pair, or the R1 hit differs markedly
# from the R2 hit, in which case both are counted as valid hits.


# print out
print(paste("Preserved all of R1 hits with:", nrow(r1_func), "hits"))
print("Number of R2 hits that were not found with paired reads in R1:")
print(nrow(r2_func_filt))
print("Making a proportion of all R2 hits:")
print(nrow(r2_func_filt)/nrow(r2_func))

# add information for each hit from which read it came from
r1_func <- cbind(r1_func, rep("R1",nrow(r1_func)))
r2_func_filt <- cbind(r2_func_filt, rep("R2", nrow(r2_func_filt)))

colnames(r1_func)[ncol(r1_func)] <- "hit_readpair"
colnames(r2_func_filt)[ncol(r2_func_filt)] <- "hit_readpair"

# concatenate all into one
colnames(r1_func) <- colnames(r2_func_filt)
all_hits <- rbind(r1_func, r2_func_filt)
print(paste("Succesfully combined R1 and R2 with", nrow(all_hits), "hits in total"))

# add RPK, RPM, RPKM for the calculation of TPM later
print("Adding count info and counts normalized in various ways (e.g. RPK)")

# each hit represent a count
all_hits$count <- rep(1, nrow(all_hits))

# calculate reads per million (RPM)
scaling_factor <- sum(all_hits$count)/10^6
all_hits$rpm <- all_hits$count/scaling_factor

# calculate read per kilobase (RPK)
all_hits$rpk <- all_hits$count/((as.numeric(all_hits$slen)*3)/1000) # reads per kilobase, length in kb is 3*prot seq length(length in amino acids) / 1000

# calculate RPKM
all_hits$rpkm <- as.numeric(all_hits$rpm)/((as.numeric(all_hits$slen)*3)/1000) 

# get RPK scaling factor, this is mainly why we wanted to estimate the total amount of prokaryotic genes in a sample. Even those with no KO annotations.
rpk_scaling_factor <- sum(as.numeric(all_hits$rpk))/10^6

print("Succesfully added various measures")

# summarize to gene level and include the KO annotation - the same genes should have the same KO annotation
print("Summarizing to gene level")
all_gene_aggregated <- aggregate(cbind(count, rpm, rpk, rpkm) ~ gene+ko, all_hits, sum)
print(paste("Successfully summarized to gene level with:", nrow(all_gene_aggregated), "genes"))

# save the combined hits
print("Saving the processed files")
setwd(args[5])

dir.create("combined_hits")
setwd("combined_hits")
fil_n <- paste(sample_namn, ".txt", sep = "")
write.table(x = all_hits, file = fil_n, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
setwd("../")

dir.create("combined_hits_gene_aggregated")
setwd("combined_hits_gene_aggregated")
write.table(x = all_gene_aggregated, file = fil_n, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
setwd("../")

# save the sample specific details, scaling factors, average gene length, etc.
dir.create("sample_details")
setwd("sample_details")
sample_details <- c(scaling_factor, rpk_scaling_factor, sum(all_hits$count), sum(all_hits$rpk), sum(all_hits$rpkm), median(all_hits$rpk))
names(sample_details) <- c("Count_PerMillion_Scaling_Factor", "RPK_PerMillion_Scaling_Factor", "Total_Counts", "Total_RPK", "Total_RPKM", "Median_RPK")
write.table(x = sample_details, file = fil_n, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

# print out session info
print("SessionInfo:")
sessionInfo()