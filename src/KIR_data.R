
## ==========================================================
##
##   KIR imputation on FinnGen array and Histogenetics data
##
## ==========================================================


## libs, options and functions
source('./src/functions.R')


## ----------------------------------------------------------
## reading and cleaning of data sets
## ----------------------------------------------------------

## KIR typing reference data (Histogenetics)

# typing results
hg.kir <- fread('./data/Histogenetics_KIR/KIR_Histogenetics_samples_all.tsv', data.table=F)
hg.kir <- gsub('+', '1', hg.kir %>% as.matrix, fixed=T)
hg.kir <- gsub('-', '0', hg.kir, fixed=T) %>% data.frame

# remove gene if zero variance
hg.kir <- hg.kir[, c(1, which(apply(hg.kir[, -1], 2, sd)>0)+1)] 

# remove duplicated rows
hg.kir <- hg.kir[!(hg.kir$SampleID %>% duplicated), ]

# sort KIR gene names alphabetically
hg.kir <- dplyr::select(hg.kir, c(colnames(hg.kir)[1], colnames(hg.kir)[-1] %>% sort))


## read genotype dosage data
fg.kir <- fread('./data/genotypes/KIR.raw', data.table=F)[, -c(2:6)] %>% filter(., FID %in% hg.kir$SampleID)
colnames(fg.kir)[-1] <- str_split_fixed(colnames(fg.kir)[-1], '_', 5) %>% data.frame %>% unite(., 'TT', 1,2,5, sep='_') %>% .$TT

## train and test subsets

# indices
train.ind <- sample(1:nrow(hg.kir), nrow(hg.kir)*0.5, F) %>% sort  
# match sample order
hg.kir.train <- hg.kir[train.ind, ]
fg.kir.train <- fg.kir[match(hg.kir.train$SampleID, fg.kir$FID), ] %>% na.omit
hg.kir.test  <- hg.kir[-train.ind, ]
fg.kir.test  <- fg.kir[match(hg.kir.test$SampleID, fg.kir$FID), ] %>% na.omit

## check overlap between train and test: should be 0
intersect(hg.kir.train$SampleID, hg.kir.test$SampleID)

## check that all samples match in geno and pheno data sets
all(hg.kir.train$SampleID==hg.kir.train$SampleID)
all(hg.kir.test$SampleID==hg.kir.test$SampleID)

## read WGS kpi KIR accuracies by Chen et al. 2020
wgs.kir <- fread('./data/1kG/KIR_WGS_accuracies.tsv', data.table=F, skip=2)
colnames(wgs.kir) <- c('KIR', 'Value')
wgs.kir$Value <- gsub('%', '', wgs.kir$Value, fixed=T) %>% as.numeric %>% "/"(100)
wgs.kir$Group <- 'WGS (data from Chen et al. 2020)'
wgs.kir$Metric <- 'Accuracy'
wgs.kir <- wgs.kir[, c('Group', 'Metric', 'KIR', 'Value')]
wgs.kir$KIR <- paste0('KIR_', wgs.kir$KIR)
wgs.kir <- filter(wgs.kir, KIR %in% kir.model.preds.acc$KIR)
