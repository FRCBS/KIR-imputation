
## ==========================================================
##
##   KIR imputation on FinnGen array and Histogenetics data
##
## ==========================================================


## libs, options and functions
source('./src/functions.R')


## ----------------------------------------------------------
## models
## ----------------------------------------------------------

## training data; fit models

kir.model.fits <- map(2:ncol(hg.kir.train), function(x) {
  fit_VS_RF(fg.kir.train, hg.kir.train[, c(1, x)], imp.thrs=5e-5) 
})
names(kir.model.fits) <- colnames(hg.kir.train)[2:ncol(hg.kir.train)]


## test data; predictions

kir.model.preds <- map(2:ncol(hg.kir.test), function(x) {
  predict_VS_RF(fg.kir.test, hg.kir.test[, c(1, x)], kir.model.fits[[x-1]]$Model) 
})
names(kir.model.preds) <- colnames(hg.kir.test)[2:ncol(hg.kir.test)]


## fit models to full data

#colnames(fg.kir)[-1] <- str_split_fixed(colnames(fg.kir)[-1], '_', 5) %>% data.frame %>% unite(., 'TT', 1,2,5, sep='_') %>% .$TT
kir.model.full <- map(2:ncol(hg.kir), function(x) {
  fit_VS_RF(fg.kir, hg.kir[, c(1, x)], imp.thrs=5e-5) 
})
names(kir.model.full) <- colnames(hg.kir)[-1]

# save models
map2(kir.model.full, names(kir.model.full), function(x, y) {
  saveRDS(x$Model, paste0("./models/", y, ".rds"))
})


## ----------------------------------------------------------
## variants used in models
## ----------------------------------------------------------

# used variants by gene

kir.model.full.gene.var <- map2(kir.model.full, kir.model.full %>% names, function(x, y) {
  data.frame(Imputation_target_gene=y, 
            ModelVariantName=names(x$Features),
            Chr=19,
            Pos=str_split_fixed(names(x$Features), '_', 3)[, 2],
            CountedAllele=str_split_fixed(names(x$Features), '_', 3)[, 3],
            Importance=x$Features
            ) %>% return()
})

# all variants into a single table
kir.model.full.vars <- c(map(kir.model.full, function(x) names(x$Features)) %>% unlist,
                         map(kir.model.fits, function(x) names(x$Features)) %>% unlist) %>% unique
kir.model.full.vars <- data.frame(Var=kir.model.full.vars, str_split_fixed(kir.model.full.vars, '_', 3)[, 1:3])
colnames(kir.model.full.vars)[2:4] <- c('Chr', 'Pos', 'Counted_allele')
saveRDS(kir.model.full.vars, './models/model_SNP_data.rds')
var.data <- readRDS('./models/model_SNP_data.rds')

# extract variant freq. means from training genotype data
kir.model.full.vars$Counted_allele_means <- dplyr::select(fg.kir, kir.model.full.vars$Var) %>% colMeans

# original variant names
fg.vars <- fread('./data/genotypes/KIR.raw', data.table=F)[, -c(1:6)] %>% colnames
fg.vars <- data.frame(Var=fg.vars, 
                      (str_split_fixed(fg.vars, '_', 5)[, c(1,2,5)] %>% data.frame %>% unite(., 'version', 1:3, sep='_')),
                      (str_split_fixed(fg.vars, '_', 5)[, c(3:4)] %>% data.frame))
fg.vars <- fg.vars[match(kir.model.full.vars$Var, fg.vars$version), ]

# is counted dosage allele minor or major
kir.model.full.vars$MinorAllele <- NA
kir.model.full.vars$MinorAllele[ifelse(kir.model.full.vars$Counted_allele_means<0.5, T, F)] <-
  kir.model.full.vars$Counted_allele[ifelse(kir.model.full.vars$Counted_allele_means<0.5, T, F)]

# add minor and major alleles
kir.model.full.vars <- data.frame(kir.model.full.vars,
                                  map(1:length(kir.model.full.vars$MinorAllele), function(i) {
                                    out <- NA
                                    if(is.na(kir.model.full.vars$MinorAllele[i])) {
                                      major <- kir.model.full.vars$Counted_allele[i]
                                      minor <- fg.vars[i, 3:4][fg.vars[i, 3:4]!=kir.model.full.vars$Counted_allele[i]]
                                      out <- c(minor, major)
                                    } else {
                                      minor <- kir.model.full.vars$MinorAllele[i]
                                      major <- fg.vars[i, 3:4][fg.vars[i, 3:4]!=kir.model.full.vars$Counted_allele[i]]
                                      out <- c(minor, major)
                                    }
                                    return(out)
                                  }) %>% do.call(rbind, .)
)
kir.model.full.vars <- dplyr::select(kir.model.full.vars, -MinorAllele)
colnames(kir.model.full.vars)[6:7] <- c('Minor', 'Major')


## ----------------------------------------------------------
## prediction of test set with 
## varying fraction of missing snps
## ----------------------------------------------------------

fg.kir.bim <- fread('./data/genotypes/KIR.bim', header=F, data.table=F)
fg.kir.raw.test <- fread('./data/genotypes/KIR.raw', data.table=F) %>% dplyr::filter(., FID %in% fg.kir.test$FID)
colnames(fg.kir.raw.test)[-c(1:6)] <- colnames(fg.kir)[-1]

# loop over fractions
kir.model.preds.missing <- map(seq(0.1, 0.99, by=0.1), function(f) {
  
  print(f)
  # repeat for a given fraction
  map(1:10, function(i) {
    
    # fraction of SNPs present
    var.fraction <- f
    # random sampling of variants
    tmp.ind <- sample(7:ncol(fg.kir.raw.test), (ncol(fg.kir.raw.test)-6)*var.fraction, F) %>% sort
    tmp.dat <- checkInputData(fg.kir.raw.test[, c(1:6, tmp.ind)], fg.kir.bim[tmp.ind-6, ], var.data)
    # prediction of test set
    kir.model.preds.tmp <- map(2:ncol(hg.kir.test), function(x) {
      predict_VS_RF(tmp.dat, hg.kir.test[, c(1, x)], kir.model.fits[[x-1]]$Model) 
    })
    names(kir.model.preds.tmp) <- colnames(hg.kir.test)[2:ncol(hg.kir.test)]
    # error calculation
    map2(kir.model.preds.tmp, names(kir.model.preds.tmp), function(x, y) { 
      out <- confusionMatrix(ifelse(x$Predicted>0.5, 1, 0) %>% factor, 
                             reference=x$Actual %>% factor)$byClass %>% data.frame
      data.frame(Group='test set', Fraction=var.fraction, Iteration=i, Metric=rownames(out), KIR=y, Value=out[, 1])
    }) %>% do.call(rbind, .)
    
  }) %>% do.call(rbind, .)
  
}) %>% do.call(rbind, .)



## ----------------------------------------------------------
## extract prediction errors
## ----------------------------------------------------------

## training data

# table for boxplots
kir.model.fits.oob <- map2(kir.model.fits, colnames(hg.kir.train)[-1], function(x, y) { 
  out <- rbind(data.frame(Genotype='present', Predicted=filter(x$OOB, Actual==1)$Predicted),
               data.frame(Genotype='absent',  Predicted=filter(x$OOB, Actual==0)$Predicted))
  data.frame(Gene=y, out)
}) %>% do.call(rbind, .)
kir.model.fits.oob$Genotype <- ifelse(kir.model.fits.oob$Genotype=='present', 1, 0)

# confusion 
kir.model.fits.conf <- map2(kir.model.fits, names(kir.model.fits), function(x, y) { 
  out <- confusionMatrix(ifelse(x$OOB$Predicted>0.5, 1, 0) %>% factor, 
                         reference=x$OOB$Actual %>% factor)$table %>% data.frame
  data.frame(KIR=y, out)
}) %>% do.call(rbind, .)
kir.model.fits.conf$Prediction <- factor(kir.model.fits.conf$Prediction, levels=c(1, 0)) # arrange row levels for plotting
# fill color for plotting
kir.model.fits.conf$Col <- rep(c('green', 'red', 'red', 'green'), kir.model.fits.conf$KIR %>% unique %>% length)

# Accuracy metrics
kir.model.fits.metrics <- map2(kir.model.fits, names(kir.model.fits), function(x, y) { 
  out <- confusionMatrix(ifelse(x$OOB$Predicted>0.5, 1, 0) %>% factor, 
                         reference=x$OOB$Actual %>% factor)$byClass %>% data.frame
  data.frame(Group='training set (OOB)', Metric=rownames(out), KIR=y, Value=out[, 1])
}) %>% do.call(rbind, .)


## test data

# table for beeswarms
kir.model.preds.tab <- map2(kir.model.preds, colnames(hg.kir.test)[-1], function(x, y) { 
  out <- rbind(data.frame(Genotype='present', Predicted=filter(x, Actual==1)$Predicted),
               data.frame(Genotype='absent',  Predicted=filter(x, Actual==0)$Predicted))
  data.frame(Gene=y, out)
}) %>% do.call(rbind, .)
kir.model.preds.tab$GenoNum <- ifelse(kir.model.preds.tab$Genotype=='present', 1, 0)
kir.model.preds.tab$Gene <- gsub('KIR_', '', kir.model.preds.tab$Gene)

# confusion 
kir.model.preds.conf <- map2(kir.model.preds, names(kir.model.preds), function(x, y) { 
  out <- confusionMatrix(ifelse(x$Predicted>0.5, 1, 0) %>% factor, 
                         reference=x$Actual %>% factor)$table %>% data.frame
  data.frame(KIR=y, out)
}) %>% do.call(rbind, .)
kir.model.preds.conf$Prediction <- factor(kir.model.preds.conf$Prediction, levels=c(1, 0)) # arrange row levels for plotting
# fill color for plotting
kir.model.preds.conf$Col <- rep(c('green', 'red', 'red', 'green'), kir.model.preds.conf$KIR %>% unique %>% length)
kir.model.preds.conf$KIR <- gsub('KIR_', '', kir.model.preds.conf$KIR)

# Accuracy metrics
kir.model.preds.metrics <- map2(kir.model.preds, names(kir.model.preds), function(x, y) { 
  out <- confusionMatrix(ifelse(x$Predicted>0.5, 1, 0) %>% factor, 
                         reference=x$Actual %>% factor)$byClass %>% data.frame
  data.frame(Group='test set', Metric=rownames(out), KIR=y, Value=out[, 1])
}) %>% do.call(rbind, .)

kir.model.preds.acc <- map2(kir.model.preds, names(kir.model.preds), function(x, y) { 
  out <- confusionMatrix(ifelse(x$Predicted>0.5, 1, 0) %>% factor, 
                         reference=x$Actual %>% factor)$overall %>% data.frame
  data.frame(Group='test set', Metric=rownames(out), KIR=y, Value=out[, 1])
}) %>% do.call(rbind, .)


## full data

# table for boxplots
kir.model.full.oob <- map2(kir.model.full, colnames(hg.kir)[-1], function(x, y) { 
  out <- rbind(data.frame(Genotype='present', Predicted=filter(x$OOB, Actual==1)$Predicted),
               data.frame(Genotype='absent',  Predicted=filter(x$OOB, Actual==0)$Predicted))
  data.frame(Gene=y, out)
}) %>% do.call(rbind, .)
kir.model.full.oob$Genotype <- ifelse(kir.model.full.oob$Genotype=='present', 1, 0)

# confusion 
kir.model.full.conf <- map2(kir.model.full, names(kir.model.full), function(x, y) { 
  out <- confusionMatrix(ifelse(x$OOB$Predicted>0.5, 1, 0) %>% factor, 
                         reference=x$OOB$Actual %>% factor)$table %>% data.frame
  data.frame(KIR=y, out)
}) %>% do.call(rbind, .)
kir.model.full.conf$Prediction <- factor(kir.model.full.conf$Prediction, levels=c(1, 0)) # arrange row levels for plotting
# fill color for plotting
kir.model.full.conf$Col <- rep(c('green', 'red', 'red', 'green'), kir.model.full.conf$KIR %>% unique %>% length)

# Accuracy metrics
kir.model.full.metrics <- map2(kir.model.full, names(kir.model.full), function(x, y) { 
  out <- confusionMatrix(ifelse(x$OOB$Predicted>0.5, 1, 0) %>% factor, 
                         reference=x$OOB$Actual %>% factor)$byClass %>% data.frame
  data.frame(Group='full data (OOB)', Metric=rownames(out), KIR=y, Value=out[, 1])
}) %>% do.call(rbind, .)


## WGS
wgs.comp <- rbind(kir.model.preds.metrics %>% filter(Metric=='Accuracy'),
                  wgs.kir)

