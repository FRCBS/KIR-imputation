
## train models for R9 SNPs

## libs, options and functions
source('./src/functions.R')


## ----------------------------------------------------------
## data
## ----------------------------------------------------------

# typing results
hg.kir <- read_xlsx('data/phenotypes/KIR_Histogenetics_samples_all.xlsx')
hg.kir <- gsub('+', '1', hg.kir %>% as.matrix, fixed=T)
hg.kir <- gsub('-', '0', hg.kir, fixed=T) %>% data.frame

# remove gene if zero variance
hg.kir <- hg.kir[, c(1, which(apply(hg.kir[, -1], 2, sd)>0)+1)] 

# remove duplicated rows
hg.kir <- hg.kir[!(hg.kir$SampleID %>% duplicated), ]

# sort KIR gene names alphabetically
hg.kir <- dplyr::select(hg.kir, c(colnames(hg.kir)[1], colnames(hg.kir)[-1] %>% sort))


## read genotype dosage data
fg.kir <- fread('data/genotypes/KIR.raw', data.table=F)[, -c(2:6)] %>% filter(., FID %in% hg.kir$SampleID)
colnames(fg.kir)[-1] <- str_split_fixed(colnames(fg.kir)[-1], '_', 5) %>% data.frame %>% unite(., 'TT', 1,2,5, sep='_') %>% .$TT

# R9 SNP content
r9 <- fread('data/genotypes/BB_finngen_R9_chr19_KIR.vcf.gz.bim', data.table = F)

# select common SNPs with R9
# returns indices of fg.kir that match to r9
find.ind <- match(r9$V2 %>% paste0('chr', .),
                  str_split_fixed(colnames(fg.kir), '_', 3)[, 1:2] %>% 
                    data.frame %>% unite(U) %>% .$U)
fg.kir.r9 <- fg.kir[, c(1, find.ind %>% na.omit)]
# check kept snps
intersect(r9$V2 %>% paste0('chr', .), 
          str_split_fixed(colnames(fg.kir.r9), '_', 3)[, 1:2] %>% data.frame %>% unite(U) %>% .$U) %>% 
  length
nrow(r9)
ncol(fg.kir)
ncol(fg.kir.r9)


## ----------------------------------------------------------
## models
## ----------------------------------------------------------

kir.model.r9 <- map(2:ncol(hg.kir), function(x) {
  fit_VS_RF(fg.kir.r9, hg.kir[, c(1, x)], imp.thrs=5e-5) 
})
names(kir.model.r9) <- colnames(hg.kir)[-1]

# save models
map2(kir.model.r9, names(kir.model.r9), function(x, y) {
  saveRDS(x$Model, paste0("./models/R9/", y, "_R9.rds"))
})

# all variants into a single table
kir.model.r9.vars <- c(map(kir.model.r9, function(x) names(x$Features)) %>% unlist) %>% unique
kir.model.r9.vars <- data.frame(Var=kir.model.r9.vars, str_split_fixed(kir.model.r9.vars, '_', 3)[, 1:3])
colnames(kir.model.r9.vars)[2:4] <- c('Chr', 'Pos', 'Counted_allele')
saveRDS(kir.model.r9.vars, './models/R9/model_SNP_data.rds')

# model variants
model.fits.vars <- map(kir.model.r9, function(x) names(x$Features)) %>% unlist %>% unique
model.fits.vars <- data.frame(Var=model.fits.vars, str_split_fixed(model.fits.vars, '_', 3)[, 1:3])
colnames(model.fits.vars)[2:4] <- c('Chr', 'Pos', 'Counted_allele')
# extract variant freq. means from training genotype data
model.fits.vars$Counted_allele_means <- dplyr::select(fg.kir.r9, model.fits.vars$Var) %>% colMeans
# save to RDS file
saveRDS(model.fits.vars, paste0(args[4], 'models/SNP_data.rds'))
# save plink allele reference file
fwrite(data.frame(unite(model.fits.vars[, 2:3], U)$U %>% gsub('chr', '', .), model.fits.vars[, 4]), 
       'models/R9/plink_allele_ref', col.names=F, row.names=F, sep='\t')



