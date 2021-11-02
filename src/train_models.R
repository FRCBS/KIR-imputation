
## libs
suppressMessages(library(ranger))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

# TO DO: produce a plink allele ref

## command line arguments:
## 1. path to genotype plink bim file (.bim)
## 2. path to genotype plink dosage file (.raw)
## 3. path to phenotype file
## 4. model output path
args <- commandArgs(trailingOnly=T)[1:4]
if(length(na.omit(args))!=4) stop("Incorrect number of input arguments: stopping.")


## ----------------------------------------------------------
## functions
## ----------------------------------------------------------

# pheno.dat <- pheno.kir[, c(1, 2)]
# geno.dat <- geno.raw
# tmp.ind <- apply(geno.dat, 2, function(x) any(is.na(x)))
# geno.dat <- geno.dat[, tmp.ind==F]
fit_VS_RF <- function(geno.dat, pheno.dat, imp.thrs=5e-5) {
  
  # input data: first column = sample ID
  colnames(geno.dat)[1] <- colnames(pheno.dat)[1] <- 'ID'
  
  # join pheno and geno data by ID
  # output: first column = KIR phenotype
  tmp.dat  <- inner_join(pheno.dat, geno.dat, by='ID')[, -1] %>% na.omit
  kir.gene <- colnames(tmp.dat)[1] # save gene name
  cat(paste0(kir.gene, '\n'))
  colnames(tmp.dat)[1] <- 'KIR'
  tmp.dat$KIR <- factor(tmp.dat$KIR)
  
  # generate class weights
  vweights  <- table(tmp.dat[, 1])/length(tmp.dat[, 1]) # pheno level freqs
  #vweights  <- vweights[c(2, 1)] # reverse level freqs
  vweights  <- 1-vweights
  
  # select important features via rf
  selected.vars <- ranger(y=tmp.dat$KIR, x=as.matrix(tmp.dat[, -1]), probability=F, importance='permutation', class.weights=vweights,
                          mtry=ncol(tmp.dat)*0.8, num.trees=2000) %>% importance 
  selected.vars <- selected.vars[selected.vars > imp.thrs] # importance threshold
  
  # select important variants from input data
  tmp.dat <- tmp.dat[, c(1, which(colnames(tmp.dat) %in% names(selected.vars)))]
  
  # train model on important variants
  model.train <- ranger(y=tmp.dat$KIR, x=as.matrix(tmp.dat[, -1]), probability=T, mtry=ncol(tmp.dat)*0.8, class.weights=vweights,
                        importance='permutation', num.trees=2000)
  
  #cat(paste0('\n'))
  
  # return trained model and selected variants 
  list(Gene=kir.gene, Model=model.train, Features=selected.vars, 
       OOB=data.frame(Predicted=model.train$predictions, Actual=tmp.dat$KIR))
  
}



## ----------------------------------------------------------
## fit models
## ----------------------------------------------------------

## read geno and pheno data


# geno.bim <- fread("./tmp/1kG_KIR.bim", data.table=F, header=F)
# geno.raw <- fread("./tmp/1kG_KIR.raw", data.table=F)[, -c(1, 3:6)] 
# pheno.kir <- fread("./test/1kG_KIR_testpheno.tsv", data.table=F) 
# rm(geno.bim); rm(geno.raw); rm(model)


cat("Reading input data...\n")

# geno data
geno.bim <- fread(args[1], data.table=F, header=F)
geno.raw <- fread(args[2], data.table=F)[, -c(1, 3:6)]

# replace illegcal characters in variant names
colnames(geno.raw) <- gsub('<', '', colnames(geno.raw), fixed=T)
colnames(geno.raw) <- gsub('>', '', colnames(geno.raw), fixed=T)

# pick the dosage allele
counted.allele <- colnames(geno.raw)[-1] %>% str_split(., '_', 10) %>% map(., function(x) x[length(x)]) %>% unlist
# new position based input variant names
input.vars <- data.frame(paste0('chr', geno.bim[, 1]), geno.bim[, 4], counted.allele) %>% 
  unite(., 'X', 1:3, sep='_') %>% .$X
colnames(geno.raw)[-1] <- input.vars

cat(paste0("Read ", length(input.vars), " variants in input genotype data.\n"))

# pheno data
pheno.kir <- fread(args[3], data.table=F)
colnames(pheno.kir)[1] <- 'IID'

# harmonize geno and pheno sample content
common.samples <- intersect(geno.raw$IID, pheno.kir$IID)
cat(paste0("Read ", length(common.samples), " samples in genotype and phenotype data.\n"))
geno.raw <- geno.raw[match(common.samples, geno.raw$IID), ]
pheno.kir <- pheno.kir[match(common.samples, pheno.kir$IID), ]

# remove missing (NA) variants
tmp.ind <- apply(geno.raw, 2, function(x) any(is.na(x)))
geno.raw <- geno.raw[, tmp.ind==F]

## fit models
cat("Fitting models:\n")
model.fits <- map(2:ncol(pheno.kir), function(x) {
  fit_VS_RF(geno.raw, pheno.kir[, c(1, x)], imp.thrs=5e-5) 
})
names(model.fits) <- colnames(pheno.kir)[2:ncol(pheno.kir)]

# save models
cat(paste0("Saving models to ", args[4], "\n"))
map2(model.fits, names(model.fits), function(x, y) {
  saveRDS(x$Model, paste0(args[4], "/model_", y, ".rds"))
}) %>% invisible()

# model variants
model.fits.vars <- map(model.fits, function(x) names(x$Features)) %>% unlist %>% unique
model.fits.vars <- data.frame(Var=model.fits.vars, str_split_fixed(model.fits.vars, '_', 3)[, 1:3])
colnames(model.fits.vars)[2:4] <- c('Chr', 'Pos', 'Counted_allele')
# extract variant freq. means from training genotype data
model.fits.vars$Counted_allele_means <- dplyr::select(geno.raw, model.fits.vars$Var) %>% colMeans
# save to RDS file
saveRDS(model.fits.vars, paste0(args[4], '/SNP_data.rds'))
# save plink allele reference file
fwrite(data.frame(unite(model.fits.vars[, 2:3], U)$U %>% gsub('chr', '', .), model.fits.vars[, 4]), 
       paste0(args[4], '/plink_allele_ref'), col.names=F, row.names=F, sep='\t')

# save OOB estimate of fitted models
errs <- map2(model.fits, names(model.fits), function(x, y) {
  pred <- apply(x$Model$predictions, 1, function(ll) colnames(x$Model$predictions)[which.max(ll)])
  err <- pred==dplyr::select(pheno.kir, all_of(y))[,1]
  err <- 1-(sum(err)/length(err))
  data.frame(Gene=y, OOB_error=err) %>% return
}) %>% do.call(rbind, .)
fwrite(errs, paste0(args[4], '/OOB_errors.tsv'), sep='\t')



