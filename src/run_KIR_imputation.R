
## libs
suppressMessages(library(ranger))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))


## command line arguments:
## 1. path to model folder
## 2. path to genotype plink dosage file (.raw)
## 3. path to genotype plink bim file (.bim)
## 4. output path
args <- commandArgs(trailingOnly=T)[1:4]
if(length(na.omit(args))!=4) stop("Incorrect number of input arguments: stopping.")


## ----------------------------------------------------------
## Functions
## ----------------------------------------------------------

# input.raw <- fread("./tmp/1kG_KIR.raw", data.table=F)
# input.bim <- fread("./tmp/1kG_KIR.bim", data.table=F, header=F)
# model.snps <- readRDS('./tmp/tmpmodels/SNP_data.rds')
# models <- map(list.files('./tmp/tmpmodels', 'model.+.rds', full.names=T), readRDS)

# SNP harmonization
# harmonizeVars <- function(snpref, bim, raw) {
#   snpref$Pos
#   map(str_split(colnames(raw)[-c(1:6)], '_'), function(x) x[length(x)])
#   bim$V4
# }

# input data integrity check and processing
checkInputData <- function(input.raw, input.bim, model.snps) {
  
  # Use allele reference file to convert to dasage format:
  # plink --bfile input_data --recode A --recode-allele ./test/plink_allele_ref
  
  # This function checks the input dosage file snps and replaces missing ones with mean genotype dosage
  
  # input genotype dosage data: separate sample names from genos
  input.samples <- input.raw$IID
  cat(paste0('Read ', length(input.samples), ' samples in the input dosage data.\n'))
  input.dat     <- input.raw[, -c(1:6)]
  # change variable names
  input.vars <- str_split(colnames(input.dat), '_') %>% do.call(rbind, .)
  input.vars <- data.frame(paste0('chr', input.bim[, 1]), input.bim[, 4], input.vars[, ncol(input.vars)]) %>% 
    unite(., 'X', 1:3, sep='_') %>% .$X
  # replace illegcal characters in variant names
  input.vars <- gsub('<', '', input.vars, fixed=T)
  input.vars <- gsub('>', '', input.vars, fixed=T)
  
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
             round(100*(length(input.ind)/length(model.snps$Var)), 2), "%) in input genotype data.\n") )
  
  # which model variants match with input variants 
  ref.ind <- match(colnames(input.dat), model.snps$Var)
  # replace into the mean value matrix the actual input values based on matching
  tmp.mat[, ref.ind] <- input.dat
  
  # replace missing values with mean
  map(1:ncol(tmp.mat), function(x) {
    nas <- is.na(tmp.mat[, x])
    if(length(nas)>0) {
      tmp.mat[nas, x] <<- model.snps[x, 'Counted_allele_means']
    }
  }) %>% invisible()
  
  # output is a data matrix with all the model variants, and value replacement with mean where no match was found
  return(data.frame(ID=input.samples, tmp.mat, stringsAsFactors=F))
  
}

# Predict new data using fitted model and seleted variants
# return class pp
impute_VS_RF <- function(geno.dat, model) {
  
  library(ranger)
  
  # variant IDs
  vars <- model$variable.importance %>% names
  
  # input data: first column = sample ID
  colnames(geno.dat)[1] <- 'ID'
  
  # select predictor variants 
  geno.dat <- geno.dat[, c(1, which(colnames(geno.dat) %in% vars))]
  
  # predict and return pp
  out <- predict(model, geno.dat[, -1])$predictions#[, 2] 
  
  # classify result
  #classes <- ifelse(out>0.5, 1, 0)
  
  # return prediction result and actual phenotype
  data.frame(sample.id=geno.dat[, 1], posterior.probability=out)
  
} 

## ----------------------------------------------------------
## impute KIRs
## ----------------------------------------------------------

# read models
models <- map(list.files(args[1], 'model.+.rds', full.names=T), readRDS)
if(length(na.omit(models))==0) stop("No models found! Stopping.")

# KIR gene names to models
model.names <- gsub('.rds', '', list.files(args[1], 'model.+.rds'), fixed=T)

# read model SNPs
model.snps <- readRDS(paste0(args[1], '/SNP_data.rds')) 
if(length(na.omit(model.snps))==0) stop("'SNP_data.rds' file not found. Stopping.")

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
  fwrite(out, paste0(args[4], '/imputed_', gsub('model_', '', y), '.tsv'), sep='\t')
}) %>% invisible

cat('Imputation done.\n')
cat(paste0('Results written to ',  args[4], '\n'))

# Rscript ./src/run_KIR_imputation.R "./results/models" "./test/simulated_ref.raw" "./test/simulated_ref.bim" "./test/output"

