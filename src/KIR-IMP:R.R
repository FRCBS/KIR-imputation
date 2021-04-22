
## ==========================================================
##
##   KIR imputation on FinnGen array and Histogenetics data
##
## ==========================================================


## libs, options and functions
source('./src/functions.R')


## ----------------------------------------------------------
## prepare data for KIR*IMP 
## ----------------------------------------------------------

# liftover to hg19
bim <- fread('./data/FG_KIR.bim', data.table=F)
write(map(bim$V4, function(x) paste0('chr19:', x, '-', x)) %>% do.call(rbind, .), './data/KIR-IMP/vars2lift.txt', ncolumns=1)
# read results
lifted <- fread('./data/KIR-IMP/hglft_genome_15a09_dc1140.bed', data.table=F, header=F)[, 1] %>% 
  str_split_fixed(., ':|-', 3) %>% .[, 1:2] %>%  data.frame
lift.fail <- fread('./data/KIR-IMP/lift_fail.txt', header=F, data.table=F)$V1 %>% .[seq(2, length(.), by=2)] %>% 
  str_split_fixed(., ':|-', 3) %>% .[, 1:2] %>%  data.frame
# snps to remove
write(bim$V2[bim$V4 %in% lift.fail$X2], './data/KIR-IMP/remove_hg19.list', ncolumns=1)
# generate new plink data files with removed unlifted snps
system('plink --bfile ./data/FG_KIR --exclude ./data/KIR-IMP/remove_hg19.list --make-bed --out ./data/KIR-IMP/FG_KIR_hg19')
# read bim
bim <- fread('./data/KIR-IMP/FG_KIR_hg19.bim', data.table=F)
# replace old positions with lifted
bim$V4 <- lifted[, 2]
fwrite(bim, './data/KIR-IMP/FG_KIR_hg19.bim', sep='\t', col.names=F)
# phase with shapeit
system(paste0('shapeit -B ./data/KIR-IMP/FG_KIR_hg19.bed ./data/KIR-IMP/FG_KIR_hg19.bim ./data/KIR-IMP/FG_KIR_hg19.fam', 
              ' --input-from 54911716 --input-to 55511245',
              ' --burn 10 --prune 10 --main 50 --window 0.5 -O ./data/KIR-IMP/FG_KIR_hg19_phased'))
# read imputed file
haps <- fread('./data/KIR-IMP/FG_KIR_hg19_phased.haps', data.table=F) # space delimited
# read KIR*IMP SNP guide
snp.guide <- fread('./data/KIR-IMP/kirimp.uk1.snp.info.csv', data.table=F)
snp.guide <- snp.guide[snp.guide$position %in% haps$V3, ]
snp.guide$isMin <- ifelse(snp.guide$allele1_frequency<0.5, T, F)

# to harmonize snp orientation according to SNP*IMP guide
map(1:nrow(snp.guide), function(i) {
  ind <- which(haps$V3==snp.guide$position[i]) # position in haps
  tmp <- haps[ind, ] # select haps row; a data frame
  is.min <- (sum(tmp[1, 6:ncol(tmp)]) / length(6:ncol(tmp)))<0.5 # check if allele1 is minor
  if(snp.guide$isMin[i]!=is.min) { # if doesn't match with guide, flip around
    haps[ind, 4:5] <<- haps[ind, 5:4]
    haps[ind, 6:ncol(haps)] <<- ifelse(haps[ind, 6:ncol(haps)]==1, 0, 1)
    if(all(haps[ind, 4:5] == snp.guide[i, c('allele0', 'allele1')])==F) { # replace also alleles if do not match
      haps[ind, 4:5] <<- snp.guide[i, c('allele0', 'allele1')]
    }
  }
  
}) %>% invisible

# write harmonized haps
fwrite(haps, './data/KIR-IMP/FG_KIR_hg19_phased_harmonized.haps', sep=' ', col.names=F)


## ----------------------------------------------------------
## check KIR*IMP results
## ----------------------------------------------------------

# results data
res <- fread('./data/KIR-IMP/imputation_results/imputations.csv', data.table=F)
res$locus <- gsub('KIR2DS4TOTAL', 'KIR2DS4', res$locus) # measure 'Total'
res$locus <- gsub('KIR3DL1ex9', 'KIR3DL1', res$locus)
res <- map(gsub('_', '', colnames(hg.kir.test)[-1]), function(x) {
 tmp <- res[res$locus==x, ] 
 col.ind <- match(gsub('KIR', 'KIR_', x), colnames(hg.kir.test))
 out <- inner_join(tmp[, c('ID_1', 'locus', 'imputedType', 'posteriorProbability')], hg.kir[, c(1, col.ind)], 
                   by=c('ID_1'='SampleID'))
 out$imputedType <- ifelse(out$imputedType>=1, 1, 0) 
 return(out)
})
names(res) <- colnames(hg.kir)[-1]

# calculate accuracy metrics
kirimp.metrics <- map2(res, names(res), function(x, y) { 
  out <- confusionMatrix(x$imputedType %>% factor, 
                         reference=x[, 5] %>% factor)$byClass %>% data.frame
  data.frame(Group='KIR*IMP', Metric=rownames(out), KIR=y, Value=out[, 1])
}) %>% do.call(rbind, .)

kirimp.acc <- map2(res, names(res), function(x, y) { 
  out <- confusionMatrix(x$imputedType %>% factor, 
                         reference=x[, 5] %>% factor)$overall %>% data.frame
  data.frame(Group='KIR*IMP', Metric=rownames(out), KIR=y, Value=out[, 1])
}) %>% do.call(rbind, .)


# SNP allele freq correlation between input and KIR*IMP reference
tmp <- fread('./data/KIR-IMP/imputation_results/alleles_snp.csv', data.table=F)
cor(tmp$allele1_frequency_input, tmp$allele1_frequency_reference)

# mean of KIR*IMP estimates for accuracy excluding haplotypes
tmp <- fread('./data/KIR-IMP/imputation_results/accuracy.csv', data.table=F)
mean(tmp$accuracy[-c(1,2)])

###########
# tmp <- haps[which(haps$V3 %in% snp.guide$position)[1:10], 1:100]
# (sum(tmp[1, 6:ncol(tmp)]) / length(6:ncol(tmp)))<0.5
# (sum(tmp[2, 6:ncol(tmp)]) / length(6:ncol(tmp)))<0.5
# (sum(tmp[3, 6:ncol(tmp)]) / length(6:ncol(tmp)))<0.5
# (sum(tmp[4, 6:ncol(tmp)]) / length(6:ncol(tmp)))<0.5


