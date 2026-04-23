args <- commandArgs(trailingOnly = TRUE)
# requires the KEGG hits for R1 and the KEGG hits for R2 as the first and second input parameter 
# the third parameter needs to be the (downloaded) mapping file (tab delimited) between KEGG genes and KO groups (only two columns, gene and ko)
# as fourth parameter if the KEGG hits should be filtered (1) or not (0) with the defined default limits.
# currently the default limits: e-value 10^-6 and bitscore 50, change if needed
# fifth parameter needs to be the output directory (needs to exist) where the subdirectories for the results are created (if they don't exist) and results are saved

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

# set read names as rownames for easier parsing
rownames(r1_func) <- r1_func$read_location
rownames(r2_func) <- r2_func$read_location

# read-in the mapping between KEGG genes and KO groups
PROKARYOTES.DAT <- read.csv(args[3], sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F)
colnames(PROKARYOTES.DAT)=c("gene", "KO")

# map the genes to KO groups
# assign KOs, get results in a list as genes may have multiple KO annotations

# R1
kos1 <- matches(r1_func$gene, PROKARYOTES.DAT$gene, all.x = TRUE, all.y = FALSE, list = TRUE, indexes = TRUE)
kos1 <- lapply(kos1, function(x) PROKARYOTES.DAT$KO[x])
names(kos1) <- r1_func$read_location

# add them to the original data frames for easier parsing
all_kos1 <- unlist(lapply(kos1, function(x) paste(x, collapse = ";")))
r1_func$ko <- all_kos1

# R2
kos2 <- matches(r2_func$gene, PROKARYOTES.DAT$gene, all.x = TRUE, all.y = FALSE, list = TRUE, indexes = TRUE)
kos2 <- lapply(kos2, function(x) PROKARYOTES.DAT$KO[x])
names(kos2) <- r2_func$read_location

all_kos2 <- unlist(lapply(kos2, function(x) paste(x, collapse = ";")))
r2_func$ko <- all_kos2

# combine the hits from paired reads
print("Combining the hits from the paired reads")

# where both mate pairs have a hit for the same read
int_reads <- intersect(r1_func$read_location, r2_func$read_location)

# pick those reads where both have a hit
r1_func_int <- r1_func[int_reads,]
r2_func_int <- r2_func[int_reads,]

# R1 reads not having a hit for the R2 pair
r1_uniq <- r1_func[-which(r1_func$read_location%in%int_reads),]

# R2 reads not having a hit for the R1 pair
r2_uniq <- r2_func[-which(r2_func$read_location%in%int_reads),]

# these are always kept

# for the reads, where both mate pairs have a hit, but the read hits for mate pairs in same order
r1_func_int <- r1_func_int[order(r1_func_int$read_location),]
r2_func_int <- r2_func_int[order(r2_func_int$read_location),]
# this works, because there is only one (best) hit for each read. So only works with one hit / read. 

# where both mate pairs have hit exactly to the same gene, always keep these
locs_gene_match <- which(r1_func_int$gene==r2_func_int$gene)

# always keep these as well, join with r1_uniq, doesn't matter if we use r1 or r2 here
if(length(locs_gene_match)>0){
  r1_uniq <- rbind(r1_uniq, r1_func_int[locs_gene_match,])
  
  # and remove from common hits for further inspection
  r1_func_int <- r1_func_int[-locs_gene_match,]
  r2_func_int <- r2_func_int[-locs_gene_match,]
}

# for the remaining reads, where both mate pairs have a hit, first see if one the hits is clearly better and save only that hit
# if there is a bitscore difference of |>15| keep only the better scoring hit
bitscore_diff <- r1_func_int$bitscore - r2_func_int$bitscore
keep_r1 <- which(bitscore_diff >= 15)
keep_r2 <- which(bitscore_diff <= (-15))

# add to r1 uniq and r2 uniq
if(length(keep_r1)>0){
  r1_uniq <- rbind(r1_uniq, r1_func_int[keep_r1,])
}

if(length(keep_r2)>0){
  r2_uniq <- rbind(r2_uniq, r2_func_int[keep_r2,])
}

# and remove from common hits for further inspection
all_better_quality_hits <- c(keep_r1, keep_r2)
if(length(all_better_quality_hits)>0){
  r1_func_int <- r1_func_int[-all_better_quality_hits,]
  r2_func_int <- r2_func_int[-all_better_quality_hits,]
}

# for the remaining common hits
# if there is a ko annotation for r1,but not r2,keep r1
# if there is a ko annotation for r2, but not r1, keep r2
kos1_int <- kos1[rownames(r1_func_int)]
kos2_int <- kos2[rownames(r2_func_int)]

# get lengths
l_r1 <- lengths(kos1_int)
l_r2 <- lengths(kos2_int)

keep_r1 <- which(l_r2==0 & l_r1>0)
keep_r2 <- which(l_r1==0 & l_r2>0)

# add to r1 uniq and r2 uniq
if(length(keep_r1)>0){
  r1_uniq <- rbind(r1_uniq, r1_func_int[keep_r1,])
}

if(length(keep_r2)>0){
  r2_uniq <- rbind(r2_uniq, r2_func_int[keep_r2,])
}

# and remove from common hits for further inspection
all_better_quality_hits <- c(keep_r1, keep_r2)
if(length(all_better_quality_hits)>0){
  r1_func_int <- r1_func_int[-all_better_quality_hits,]
  r2_func_int <- r2_func_int[-all_better_quality_hits,]
}

# for those with ko hits in both, and there is some overlap, keep r1 hits but
# keep only the overlapping ko(s) (most conservative), discard the rest ko annotations
# for those with ko hits both, but no overlap, drop the hit entirely as unreliable
kos1_int <- kos1[rownames(r1_func_int)]
kos2_int <- kos2[rownames(r2_func_int)]

# get lengths
l_r1 <- lengths(kos1_int)
l_r2 <- lengths(kos2_int)

kos_in_both <- which(l_r1>0 & l_r2 >0)
if(length(kos_in_both)>0){
  kos1_explore <- kos1_int[kos_in_both]
  kos2_explore <- kos2_int[kos_in_both]
  new_kos <- list()
  
  for(ko in 1:length(kos1_explore)){
    inter <- intersect(kos1_explore[[ko]], kos2_explore[[ko]])
    if (length(inter) > 0) {
      kos_keep <- inter
      new_kos[[ko]] <- kos_keep
    } else { # no intersect
      new_kos[[ko]] <- NA
    }
  }
  names(new_kos) <- names(kos1_explore)
  
  # drop those hits with no intersection of kos for mate pair hits where both have ko annotations
  drop_hits <- which(is.na(new_kos))
  keep_hits <- which(!is.na(new_kos))
  
  if(length(keep_hits)>0){
    # add to r1 uniq but change ko annotation to intersection
    new_kos <- new_kos[keep_hits]
    r1_temp <- r1_func_int[kos_in_both[keep_hits],] # keep hits where both have ko.
    new_kos <- new_kos[rownames(r1_temp)]
    
    # parse
    new_kos <- unlist(lapply(new_kos, function(x) paste(x, collapse = ";")))
    r1_temp$ko <- new_kos
    
    # add to r1 uniq
    r1_uniq <- rbind(r1_uniq, r1_temp)
  }
  
  # and remove these from common hits for further inspection
  r1_func_int <- r1_func_int[-kos_in_both,]
  r2_func_int <- r2_func_int[-kos_in_both,]
}

# now we should only be left with KEGG gene hits with no ko annotation for r1 and r2
# for those hits with ko annotations in neither r1 or r2, compare the functional annotations.
# if the functional annotation is less than 0.8 similar, drop that hit entirely as unreliable. 
# for similar ones, keep r1 hits. These are only used for estimating the RPK sample sums.

# extract the reported function from stitle for these hits for r1 and r2
r1_function <- unlist(lapply(X = r1_func_int$stitle, FUN =  function(x) {
  y <- strsplit(x = x[1], split = " ")[[1]]
  paste(y[3:length(y)], collapse = " ")
}))
r2_function <- unlist(lapply(X = r2_func_int$stitle, FUN =  function(x) {
  y <- strsplit(x = x[1], split = " ")[[1]]
  paste(y[3:length(y)], collapse = " ")
}))

# put everything to lowercase
r1_function <- tolower(r1_function)
r2_function <- tolower(r2_function)

# calculate functional similarity for these hits
func_similarity <- stringsim(a = r1_function, b = r2_function)
# osa similarity, the default method is used.
# The Levenshtein distance (method='lv') counts the number of deletions, insertions and substitutions necessary to turn b into a. 
# This method is equivalent to R's native adist function.
# The Optimal String Alignment distance (method='osa') is like the Levenshtein distance but also allows transposition of adjacent 
# characters. Here, each substring may be edited only once. (For example, a character cannot be transposed twice to move it forward in the string).

drop_hits <- which(func_similarity<0.7)
# if func similarity if < 0.7 the functional descriptions are very different already
# remove all those, although some "proper hits" to close homologs for r1 and r2 might be included in this
# length difference could be explored and used somehow here?

# drop these hits which are very unsimilar in their stitle
if(length(drop_hits)>0){
  r1_func_int <- r1_func_int[-drop_hits,]
  
  # and add the remaining to r1 uniq
  r1_uniq <- rbind(r1_uniq, r1_func_int)
}

# now in r1_uniq and in r2_uniq, we have all the hits we want to keep

# concatenate all into one
all_hits <- rbind(r1_uniq, r2_uniq)
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