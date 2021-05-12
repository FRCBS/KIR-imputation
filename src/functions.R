
## ==========================================================
##
##   KIR imputation on FinnGen array and Histogenetics data
##
## ==========================================================


## libs
library(data.table)
library(tidyverse)
library(ggpubr)
library(ranger)
library(ROCR)
library(pROC)
library(caret)
library(gridExtra)
library(ggbeeswarm)
library(genbankr)
library(ggsci)
library(cowplot)

## opts
options(stringsAsFactors=F)



## ----------------------------------------------------------
## functions
## ----------------------------------------------------------

## select features and train RF model on important variants

fit_VS_RF <- function(geno.dat, pheno.dat, imp.thrs=5e-5) {
  
  # input data: first column = sample ID
  colnames(geno.dat)[1] <- colnames(pheno.dat)[1] <- 'ID'
  
  # join pheno and geno data by ID
  # output: first column = KIR phenotype
  tmp.dat  <- inner_join(pheno.dat, geno.dat, by='ID')[, -1] %>% na.omit
  kir.gene <- colnames(tmp.dat)[1]; print(kir.gene) # save gene name
  colnames(tmp.dat)[1] <- 'KIR'
  
  # generate class weights
  vweights  <- table(tmp.dat[, 1])/length(tmp.dat[, 1]) %>% unname # level freqs
  vweights  <- vweights[c(2, 1)] # reverse level freqs
  
  # select important features via rf
  selected.vars <- ranger(factor(KIR) ~., data=tmp.dat, probability=F, importance='permutation', class.weights=vweights,
                          mtry=ncol(tmp.dat)*0.8, num.trees=2000) %>% importance 
  selected.vars <- selected.vars[selected.vars > imp.thrs] # importance threshold
  
  # select important variants from input data
  tmp.dat <- tmp.dat[, c(1, which(colnames(tmp.dat) %in% names(selected.vars)))]
  
  # train model on important variants
  model.train <- ranger(factor(KIR) ~., data=tmp.dat, probability=T, mtry=ncol(tmp.dat)*0.8, class.weights=vweights,
                        importance='permutation', num.trees=2000)
  
  # return trained model and selected variants and other stuff
  list(Gene=kir.gene, Model=model.train, Features=selected.vars, 
       OOB=data.frame(Predicted=model.train$predictions[, 2], Actual=tmp.dat$KIR))
  
}


## predict new data using fitted model and seleted variants;
## return also true values

predict_VS_RF <- function(geno.dat, pheno.dat, model) {
  
  # variant IDs
  vars <- model$variable.importance %>% names
  
  # input data: first column = sample ID
  colnames(geno.dat)[1] <- colnames(pheno.dat)[1] <- 'ID'
  
  # select predictor variants 
  geno.dat <- geno.dat[, c(1, which(colnames(geno.dat) %in% vars))]
  
  # join pheno and geno data by ID
  # output: first column = KIR phenotype
  tmp.dat   <- full_join(pheno.dat, geno.dat, by='ID')[, -1] %>% na.omit
  
  model.test <- predict(model, tmp.dat[, -1])$predictions[, 2] #%>% rowMeans # PP
  
  # return prediction result and actual phenotype
  data.frame(Predicted=model.test, Actual=tmp.dat[, 1])
  
} 


## predict new data using fitted model and seleted variants
## return class-wise pp

impute_VS_RF <- function(geno.dat, model) {
  
  library(ranger)
  
  # variant IDs
  vars <- model$variable.importance %>% names
  
  # input data: first column = sample ID
  colnames(geno.dat)[1] <- 'ID'
  
  # select predictor variants 
  geno.dat <- geno.dat[, c(1, which(colnames(geno.dat) %in% vars))]
  
  # predict and return pp
  out <- predict(model, geno.dat[, -1])$predictions[, 2] #%>% rowMeans # PP
  
  # classify result
  classes <- ifelse(out>0.5, 1, 0)
  
  # pp for class=0
  out <- ifelse(classes==0, 1-out, out)
  
  # return prediction result and actual phenotype
  data.frame(sample.id=geno.dat[, 1], gene=classes, posterior.probability=out)
  
} 

# extract unoque variants from model list
# returns data frame with position and counted allele in separate columns
extractModelVars <- function(x) {
  out <- map(tmp, function(x) names(x$variable.importance)) %>% unlist %>% unique 
  out <- data.frame(Var=out, str_split_fixed(out, '_', 3)[, 1:3])
  out$X2 <- as.integer(out$X2)
 
  return(out)
}


## input data integrity check and processing
checkInputData <- function(input.raw, input.bim, model.snps) {
  
  # Use allele reference file to convert to dasage format:
  # plink --bfile input_data --recode A --recode-allele ./test/plink_allele_ref
  
  # This function checks the input dosage file snps and replaces missing ones with mean genotype dosage
  
  # input genotype dosage data: separate sample names from genos
  input.samples <- input.raw$IID
  input.dat     <- input.raw[, -c(1:6)]
  # change variable names
  input.vars <- str_split(colnames(input.dat), '_') %>% do.call(rbind, .)
  input.vars <- data.frame('chr19', input.bim[, 4], input.vars[, ncol(input.vars)]) %>% unite(., 'X', 1:3, sep='_') %>% .$X
  
  colnames(input.dat) <- input.vars
  
  # generate a matrix with mean genotype values that can be replaced by real values if they exist
  tmp.mat <- sapply(model.snps$Counted_allele_means, function(x) rep(x, nrow(input.dat))) %>% 
    data.frame(stringsAsFactors=F)
  # col. names are variant names
  colnames(tmp.mat) <- model.snps$Var
  
  # which input variants match with model variants
  input.ind <- match(model.snps$Var, colnames(input.dat)) %>% na.omit
  # keep only matching variants
  input.dat   <- input.dat[, input.ind]
  cat(paste0('Found ', length(input.ind), " shared variants (", 
             round(100*(length(input.ind)/length(model.snps$Var)), 2), "%)") )
  
  # which model variants match with input variants 
  ref.ind <- match(colnames(input.dat), model.snps$Var)
  # replace into the mean value matrix the actual input values based on matching
  tmp.mat[, ref.ind] <- input.dat
  
  # output is a data matrix with all the model variants, and value replacement with mean where no match was found
  return(data.frame(ID=input.samples, tmp.mat, stringsAsFactors=F))
  
}




