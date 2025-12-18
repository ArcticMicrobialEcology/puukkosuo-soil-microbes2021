# a small script to perform the Multiomics factor analysis (MOFA) for the Puukkosuo microbial data
# this script is run on a personal computer instead of CSC Puhti 

# requires a working python 3.11 installation on the computer and an installation of MOFA
# also requires certain R packages, such as MOFA2 and reticulate (to connect R to python), see below for the specific packages
# the analysis is run on R version 4.4.1
# mofapy version is 0.7.1

# define the version of python used for MOFA

library(reticulate)

# python directory, e.g. "C:/Program Files/Python311"
use_python(python = "")

# load needed libraries
library(MOFA2)
library(ggplot2)
library(graphics)
library(grid)
library(gridExtra)
library(randomForest)
library(caret)
library(vegan)
library(rgl)
library(lmerTest)
library(phyloseq)

# define some custom functions to be used
plotMOFA2Dord <- function(x_to_plot, y_to_plot, factor_values, metadata, colors=c("tan1", "forestgreen"), legend_graz_loc="topleft", legend_veg_loc="topright", legend_treat_loc="bottomleft"){
  
  # pick the relevant metadata
  meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
  colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")
  
  # set levels for metadata
  meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
  meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
  meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))
  
  # plot the base plot
  par(mar=c(5.1, 5.1, 4.1, 2.1))
  plot(factor_values[, x_to_plot], factor_values[, y_to_plot],
       xlab = colnames(factor_values)[x_to_plot], 
       ylab = colnames(factor_values)[y_to_plot],
       main = "Multi-omics factor analysis (MOFA)",
       bty = "n",
       cex.axis = 1.5,
       cex.lab = 1.5,
       type = "n") # create an empty plot
  
  # define colors - grazing
  cols <- rep(colors[1], nrow(meta_use))
  cols[which(meta_use$grazing=="ungrazed")] <- colors[2]
  
  # define point shapes - vegetation
  pchs <- rep(19, nrow(meta_use))
  pchs[which(meta_use$veg_clusters=="c_cho")] <- 15
  pchs[which(meta_use$veg_clusters=="c_ros")] <- 17
  
  # define point sizes - snow treatment
  cexes <- rep(2, nrow(meta_use))
  cexes[which(meta_use$treatment=="-S")] <- 1
  cexes[which(meta_use$treatment=="+S")] <- 3
  
  # add points for the samples to the plot
  points(factor_values[,x_to_plot], factor_values[,y_to_plot], col = cols, pch = pchs, cex=cexes) 
  
  # add legends
  legend(legend_graz_loc, legend = c("grazed", "ungrazed"), fill = colors, cex = 1.5)
  legend(legend_veg_loc, legend = c("t_ces", "c_cho", "c_ros") , col="black", pch = c(19,15,17), cex=1.5)
  legend(legend_treat_loc, legend = c("-S", "CTL", "+S"), col="black", pch = 19, pt.cex=c(1,2,3))
}

# load the needed datas, local copies downloaded from CSC Puhti 
# put the datas into a list for MOFA
datas_list <- list()

# kegg KO MG data - processed filtered by the read_based processing scripts - need to be downloaded in the subdirectory (of the working directory) kegg_mg 
load("kegg_mg/Matrices_For_Downstream.RData")
# log2 transform the data for MOFA 
data_for_mofa <- log2(ko_tpm_data+1)
# add data source into rownames, since there might same names in MG and MT (e.g. same KO groups or taxons), which create problems with MOFA
rownames(data_for_mofa) <- paste(rownames(data_for_mofa), "_MG", sep = "")
datas_list[[1]] <- data_for_mofa

# phyloflash MG taxonomy data - processed filtered by the read_based processing scripts - need to be downloaded in the subdirectory (of the working directory) tax_mg 
load("tax_mg/Otu_Tax_Tables_Parsed_Filtered.RData")
# use CLR transformed NTU/OTU data for MOFA
data_for_mofa <- otu_table_clr_filt
# add data source into rownames, since there might same names in MG and MT (e.g. same KO groups or taxons), which create problems with MOFA
rownames(data_for_mofa) <- paste(rownames(data_for_mofa), "_MG", sep = "")
datas_list[[2]] <- data_for_mofa

# kegg KO MT data - processed filtered by the read_based processing scripts - need to be downloaded in the subdirectory (of the working directory) kegg_mt 
load("kegg_mt/Matrices_For_Downstream.RData")
# log2 transform the data for MOFA 
data_for_mofa <- log2(ko_tpm_data+1)
# add data source into rownames, since there might same names in MG and MT (e.g. same KO groups or taxons), which create problems with MOFA
rownames(data_for_mofa) <- paste(rownames(data_for_mofa), "_MT", sep = "")
datas_list[[3]] <- data_for_mofa

# phyloflash MT taxonomy data - processed filtered by the read_based processing scripts - need to be downloaded in the subdirectory (of the working directory) tax_mt
load("tax_mt/Otu_Tax_Tables_Parsed_Filtered.RData")
# use CLR transformed NTU/OTU data for MOFA
data_for_mofa <- otu_table_clr_filt
# add data source into rownames, since there might same names in MG and MT (e.g. same KO groups or taxons), which create problems with MOFA
rownames(data_for_mofa) <- paste(rownames(data_for_mofa), "_MT", sep = "")
datas_list[[4]] <- data_for_mofa

# check how many features each dataset has
n_features_datas <- unlist(lapply(datas_list, nrow))
n_features_datas

# try to balance the number of features in the datas as recommended in the MOFA manual.
# However, cannot fully do that, as tax mg is much less features than the other datas.

# include only 50% of the features for the KEGG KO MG and phyloflash NTU/OTU MT datasets to balance
n_included <- c(0.5,1,1,0.5)

# choose (highly) variable features for the datasets where not all features are included, as advised in the MOFA documentation
for(i in 1:length(datas_list)){
  dat <- datas_list[[i]]
  
  # use coefficient of variation to determine the highly varying features
  r_vars <- apply(dat, 1, function(x) sd(x)/mean(x))
  r_vars <- r_vars[order(r_vars, decreasing = T)]
  
  # take the top %
  included_feat <- nrow(dat) * n_included[i]
  r_vars <- r_vars[1:included_feat]
  dat <- dat[names(r_vars),]
  datas_list[[i]] <- dat
}

# check the number of features in the datas again
n_features_datas <- unlist(lapply(datas_list, nrow))
n_features_datas
# much more balanced now

# start running MOFA 
# transform the data.frames in the datas_list to matrices for MOFA (required)
for(i in 1:length(datas_list)){
  datas_list[[i]] <- as.matrix(datas_list[[i]])
}

# name the datas
names(datas_list) <- c("kegg_mg", "tax_mg", "kegg_mt", "tax_mt")

# create the MOFA object
mofaobj <- create_mofa(datas_list)

# explore the MOFA object
plot_data_overview(mofaobj)

# train MOFA
# default data options
data_opts <- get_default_data_options(mofaobj)
data_opts$scale_views <- TRUE # scale the datas as recommended by the MOFA documentation, the datas have pretty different sds

# default model options
model_opts <- get_default_model_options(mofaobj) 

# training options
train_opts <- get_default_training_options(mofaobj)

# use slow mode for best results
train_opts$convergence_mode <- "slow" 

# define random seed
train_opts$seed <- 1

# train in silent mode
train_opts$verbose <- FALSE

# prepare the MOFA object to be trained
trained_mofa_slow <- prepare_mofa(
  object = mofaobj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# MOFA says:
# In prepare_mofa(object = mofaobj, data_options = data_opts, model_options = model_opts,  :
# The total number of samples is very small for learning 15 factors.  
# Try to reduce the number of factors to obtain meaningful results. It should not exceed ~9.

# use 8 factors then
model_opts$num_factors <- 8

# prepare the MOFA object to be trained
trained_mofa_slow <- prepare_mofa(
  object = mofaobj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# everything should be good to go
# train MOFA and save the model
outfile  <- file.path(getwd(),"mofa_hvf_slow.hdf5")
trained_mofa_slow <- run_mofa(trained_mofa_slow, outfile, use_basilisk = FALSE)
saveRDS(trained_mofa_slow,"mofa_hvf_slow.rds")
# factor 1 seems to be correlating with the total number of expressed features for at least one of the datas.

# add metadata
metadata$sample <- rownames(metadata)
samples_metadata(trained_mofa_slow) <- metadata

# plot explained variance
plot_variance_explained(trained_mofa_slow, max_r2=50, las=2)
plot_variance_explained(trained_mofa_slow, max_r2=50, plot_total = TRUE)

# correlation between factors
plot_factor_cor(trained_mofa_slow)

# correlate with metadata - not stricly useful since factors which are converted to numeric but gives an idea
correlate_factors_with_covariates(trained_mofa_slow, 
                                  covariates = c("Grazing","Veg_clusters","Treatment"), 
                                  plot="log_pval")


correlate_factors_with_covariates(trained_mofa_slow, 
                                  covariates = c("Grazing","Veg_clusters","Treatment"), 
                                  plot="r")

# extract the factor values from MOFA
factor_values <- get_factors(object = trained_mofa_slow)
factor_values <- data.frame(factor_values$group1, stringsAsFactors = F, check.names = F)
factor_values <- factor_values[rownames(metadata),]

# look more specifically how the different factors are related to the exclusion treatment - to decide which factors to visualize
# simple linear regression
# and linear mixed models
# and non parametric tests

lm_models <- list()
lmm_models <- list()

grazing_p_vals <- data.frame(matrix(nrow = 3, ncol = ncol(factor_values)))
rownames(grazing_p_vals) <- c("lm", "lmm", "wilcoxon")
colnames(grazing_p_vals) <- colnames(factor_values)

# go through the factors one by one
for(i in 1:ncol(factor_values)){
  
  # gather the relevant factor values and metadata
  temp_data <- data.frame(factor_values[,i], metadata$Grazing, metadata$Treatment, metadata$Veg_clusters)
  colnames(temp_data) <- c("fac_val", "grazing", "treatment", "veg_clusters")
  rownames(temp_data) <- rownames(metadata)
  
  # set levels 
  temp_data$grazing <- factor(temp_data$grazing, levels = c("ungrazed", "grazed"))
  temp_data$treatment <- factor(temp_data$treatment, levels = c("CTL", "+S", "-S"))
  temp_data$veg_clusters <- factor(temp_data$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))
  
  # regression for factor values against metadata variables
  lm_fac_vals <- lm(fac_val ~ ., data = temp_data)
  lm_models[[i]] <- lm_fac_vals
  
  # get linear regression significance values
  sum_lm_fac_vals <- summary(lm_fac_vals)
  grazing_p_vals[1,i] <- sum_lm_fac_vals$coefficients[2,4]
  
  # linear regression using mixed modelling and setting vegetation cluster as random effect
  lmm_fac_vals <- lmer(fac_val ~ grazing + treatment + (1|veg_clusters), data = temp_data)
  lmm_models[[i]] <- lmm_fac_vals
  
  # get linear mixed effects regression significance values
  sum_lmm_fac_vals <- summary(lmm_fac_vals)
  grazing_p_vals[2,i] <- sum_lmm_fac_vals$coefficients[2,5]
  
  # perform non-parametric test
  graz_wilcox_fac_vals <- wilcox.test(x = temp_data$fac_val[which(temp_data$grazing=="grazed")], y = temp_data$fac_val[which(temp_data$grazing=="ungrazed")], alternative = "two.sided")
  
  # gather significance values
  grazing_p_vals[3,i] <- graz_wilcox_fac_vals$p.value
}

# factors 3 and 6 seem to be most (and really the only ones) associated with the exclusion treatment (grazing)

# plot these factors
{
  pdf(file = "factor3_factor6_grazing.pdf", width = 10, height = 10, onefile = T)
  plotMOFA2Dord(x_to_plot = 3, y_to_plot = 6, factor_values = factor_values, metadata = metadata, legend_graz_loc = "topright", legend_veg_loc = "bottom", legend_treat_loc = "bottomright")
  dev.off()
}

# predict the exclusion treatment (grazing here), snow treatment and vegetation cluster status for each sample using MOFA and random forest

# pick the relevant metadata
meta_use <- metadata[,c("Grazing","Treatment","Veg_clusters")]
colnames(meta_use) <- c("grazing", "treatment", "veg_clusters")

# set levels for metadata
meta_use$grazing <- factor(meta_use$grazing, levels = c("ungrazed", "grazed"))
meta_use$treatment <- factor(meta_use$treatment, levels = c("CTL", "+S", "-S"))
meta_use$veg_clusters <- factor(meta_use$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# an empty data frame where the predicted values will be saved
meta_predicted <- data.frame(matrix(nrow = nrow(meta_use), ncol = ncol(meta_use)))
rownames(meta_predicted) <- rownames(meta_use)
colnames(meta_predicted) <- colnames(meta_use)

for(i in 1:nrow(meta_use)){
  for(j in 1:ncol(meta_use)){
    
    # training data 
    temp_data <- data.frame(factor_values, meta_use[,j], check.names = F, stringsAsFactors = F)
    colnames(temp_data)[ncol(temp_data)] <- colnames(meta_use)[j]
    train_data <- temp_data[-i,]
    
    # train the model
    set.seed(i)
    if(j==1){
      model_rf <- randomForest(grazing ~ ., data=train_data, ntree=1000)
    }else if(j==2){
      model_rf <- randomForest(treatment ~ ., data=train_data, ntree=1000)
    }else{
      model_rf <- randomForest(veg_clusters ~ ., data=train_data, ntree=1000)
    }
    
    # test data
    test_data <- temp_data[i,1:(ncol(temp_data)-1)]
    
    # predict
    meta_predicted[i,j] <- as.character(stats::predict(model_rf, newdata = test_data))
  }
}

# turn into factors
for(i in 1:ncol(meta_predicted)){meta_predicted[,i] <- as.factor(meta_predicted[,i])}

# set levels matching the original metadata
meta_predicted$grazing <- factor(meta_predicted$grazing, levels = c("ungrazed", "grazed"))
meta_predicted$treatment <- factor(meta_predicted$treatment, levels = c("CTL", "+S", "-S"))
meta_predicted$veg_clusters <- factor(meta_predicted$veg_clusters, levels = c("t_ces", "c_cho", "c_ros"))

# compare - get confusion matrices for each factor, exclusion treatment (grazing), snow treatment and vegetations cluster
conf_mat_grazing <- confusionMatrix(data=meta_predicted$grazing, reference = meta_use$grazing)
conf_mat_treatment <- confusionMatrix(data=meta_predicted$treatment, reference = meta_use$treatment)
conf_mat_veg_clusters <- confusionMatrix(data=meta_predicted$veg_clusters, reference = meta_use$veg_clusters)

# gather balanced accuracies
bal_accuracies <- c(conf_mat_grazing$byClass["Balanced Accuracy"], as.vector(conf_mat_treatment$byClass[,"Balanced Accuracy"]),
                    as.vector(conf_mat_veg_clusters$byClass[,"Balanced Accuracy"]))
names(bal_accuracies) <- c("Grazing", "CTL", "+S", "-S", "T.Ces", "C.Cho", "C.Ros")

# plot the balanced accuracy
{
  pdf("balanced_accuracies.pdf", width = 11.7, height = 5)
  barplot(bal_accuracies, ylab = "Balanced Accuracy", col = "firebrick2", ylim = c(0,0.8))
  dev.off()
}

# save everything
ses_info <- sessionInfo()
save.image("mofa_done.RData")

# print out session info
print("SessionInfo:")
sessionInfo()

# postprocess in inkscape like all figures