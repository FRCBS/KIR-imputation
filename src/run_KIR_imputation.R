
## libs
suppressMessages(library(ranger))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))


## command line arguments:
## 1. path to model folder
## 2. path to genotype plink dosage file
## 3. path to genotype plink bim file
## 4. output path
args <- commandArgs(trailingOnly=T)[1:4]

if(length(na.omit(args))!=4) stop("Incorrect number of input arguments: stopping.")

## Functions

# input data integrity check and processing
checkInputData <- function(input.raw, input.bim, model.snps) {
  
  # Use allele reference file to convert to dasage format:
  # plink --bfile input_data --recode A --recode-allele ./test/plink_allele_ref
  
  # This function checks the input dosage file snps and replaces missing ones with mean genotype dosage
  
  # input genotype dosage data: separate sample names from genos
  input.samples <- input.raw$IID
  cat(paste0('Read ', length(input.samples), ' samples in the input dosage data...\n'))
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
  cat(paste0('Found ', length(input.ind), " model variants (", 
             round(100*(length(input.ind)/length(model.snps$Var)), 2), "%)...\n") )
  
  # which model variants match with input variants 
  ref.ind <- match(colnames(input.dat), model.snps$Var)
  # replace into the mean value matrix the actual input values based on matching
  tmp.mat[, ref.ind] <- input.dat
  
  # output is a data matrix with all the model variants, and value replacement with mean where no match was found
  return(data.frame(ID=input.samples, tmp.mat, stringsAsFactors=F))
  
}

# Predict new data using fitted model and seleted variants
# return class-wise pp
impute_VS_RF <- function(geno.dat, model) {
  
  library(ranger)
  
  # variant IDs
  vars <- model$variable.importance %>% names
  
  # input data: first column = sample ID
  colnames(geno.dat)[1] <- 'ID'
  
  # select predictor variants 
  geno.dat <- geno.dat[, c(1, which(colnames(geno.dat) %in% vars))]
  
  # predict and return pp
  out <- predict(model, geno.dat[, -1])$predictions[, 2] 
  
  # classify result
  classes <- ifelse(out>0.5, 1, 0)
  
  # pp for class=0
  #out <- ifelse(classes==0, 1-out, out)
  
  # return prediction result and actual phenotype
  data.frame(sample.id=geno.dat[, 1], gene=classes, posterior.probability=out)
  
} 


## impute KIRs

# read models
models <- map(list.files(args[1], 'KIR.+.rds', full.names=T), readRDS)
if(length(na.omit(models))==0) stop("No models found! Stopping.")

# KIR gene names to models
model.names <- gsub('.rds', '', list.files(args[1], 'KIR.+.rds'), fixed=T)

# read model SNPs
model.snps <- readRDS(paste0(args[1], '/model_SNP_data.rds')) 
if(length(na.omit(model.snps))==0) stop("'model_SNP_data.rds' file not found. Stopping.")

# read dosages
if(file.exists(args[2])) {
  tmp.dosages <- fread(args[2], data.table=F)
} else stop("Plink dosage file not found. Stopping.")

# read bim
if(file.exists(args[3])) {
  tmp.bim <- fread(args[3], data.table=F)
} else stop("Plink bim file not found. Stopping.")

# format input data
tmp.dosages <- checkInputData(tmp.dosages, tmp.bim, model.snps)

# predict
map2(models, model.names, function(x, y) {
  out <- impute_VS_RF(tmp.dosages, model=x)
  fwrite(out, paste0(args[4], '/imputed_', y, '.tsv'), sep='\t')
}) %>% invisible

cat('Imputation done.\n')
cat(paste0('Results written to ',  args[4], '\n'))

# Rscript ./src/run_KIR_imputation.R "./results/models" "./test/simulated_ref.raw" "./test/simulated_ref.bim" "./test/output"

