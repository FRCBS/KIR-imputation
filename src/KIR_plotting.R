## ==========================================================
##
##   KIR imputation on FinnGen array and Histogenetics data
##
## ==========================================================


## libs, options and functions
source('./src/functions.R')


## ----------------------------------------------------------
## gather all accuracy metrics data for plotting
## ----------------------------------------------------------

all.metrics        <- rbind(kir.model.fits.metrics, kir.model.preds.metrics, kir.model.full.metrics, 
                            wgs.comp, kirimp.metrics, kirimp.acc)
all.metrics$KIR    <- gsub('KIR_', '', all.metrics$KIR) 
all.metrics$KIR    <- factor(all.metrics$KIR, levels=unique(all.metrics$KIR))
all.metrics$Metric <- factor(all.metrics$Metric, levels=unique(all.metrics$Metric))


## ----------------------------------------------------------
## result plots
## ----------------------------------------------------------

## metrics

p.all.metrics <- ggplot(filter(all.metrics, 
                               Metric %in% c('Pos Pred Value', 'Neg Pred Value', 'Balanced Accuracy'),
                               Group %in% c('full data (OOB)', 'test set', 'training set (OOB)', 'KIR*IMP')), 
                        aes(KIR, Value, color=Group, group=Group)) +
  geom_point() + geom_line() + 
  scale_color_manual(values=c('#619CFF', 'grey70', '#F8766D' , '#00BA38')) +
  facet_wrap(~ Metric, nrow=4, ncol=1, scales='free_y') +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position='top', panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA, size=0.3))
p.all.metrics


## missing

kir.model.preds.missing.stats <- filter(kir.model.preds.missing, Metric %in% c('Pos Pred Value', 'Neg Pred Value', 'Balanced Accuracy')) %>% 
  group_by(., Fraction, Metric) %>% 
  summarise(., Mean=mean(Value, na.rm=T), SD=sd(Value, na.rm=T), CI_l=quantile(Value, 0.025, na.rm=T), CI_h=quantile(Value, 0.975, na.rm=T))
kir.model.preds.missing.stats$YMax <- kir.model.preds.missing.stats$Mean+kir.model.preds.missing.stats$SD
kir.model.preds.missing.stats$YMax[kir.model.preds.missing.stats$YMax>1] <- 1
kir.model.preds.missing.stats$YMin <- kir.model.preds.missing.stats$Mean-kir.model.preds.missing.stats$SD

p.kir.model.preds.missing.stats <- ggplot(kir.model.preds.missing.stats, aes(Fraction, Mean, ymax=YMax, ymin=YMin, group=1)) +
  geom_line(size=0.4, color='#F8766D') +
  geom_ribbon(aes(ymin=YMin, ymax=YMax), alpha=0.12) +
  xlab('Fraction of SNPs present') +
  ylab('Mean value ± S.D') +
  ylim(c(0.4, 1)) +
  facet_wrap(~Metric) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position='top', axis.text.x=element_text(size=6), panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA, size=0.3))


## WGS/kpi comparison

p.wgs.comp <- ggplot(filter(all.metrics, 
                            Group %in% c('test set', 'WGS (data from Chen et al. 2020)'),
                            Metric=='Accuracy'), 
                     aes(KIR, Value, group=Group, color=Group)) +
  geom_point(position=position_dodge(width=0.2)) + geom_line() +
  scale_color_manual(values=c('#F8766D', ggsci::pal_npg()(5)[4])) +
  scale_y_continuous(limits=c(0.9, 1.01), expand=c(0, 0)) +
  ylab('Overall accuracy') +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position='top', panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA, size=0.3), strip.background=element_rect(fill='grey92', color='grey92'))


## test set PP distributions and confusion tables

# post probs
p.test.bee <- ggplot(kir.model.preds.tab, aes(GenoNum, Predicted)) +
  geom_quasirandom(dodge.width=0.8, shape=21) +
  geom_hline(yintercept=0.5, color="grey", linetype="dashed", size=0.4) +
  ylab('Posterior probability') + xlab('Genotype') +
  facet_wrap(~ Gene) +
  theme_minimal() +
  theme(panel.border=element_rect(fill=NA, size=0.3), panel.grid=element_blank())

# confusion tables
p.test.conf <- ggplot(kir.model.preds.conf, aes(Reference, Prediction)) +
  geom_tile(aes(fill=Col), alpha=0.2, colour='black') +
  scale_fill_identity() + 
  geom_text(aes(label=Freq)) +
  xlab('Genotype') +
  facet_wrap(~KIR) +
  theme_minimal() +
  theme(panel.grid=element_blank())



## error rates by PP threshold

balanced.acc.pp <- map(seq(0.05, 0.4, by=0.01), function(tt) {
  pp.thres <- tt
  map2(kir.model.preds, kir.model.preds %>% names, function(x, y) {
    out <- filter(x, Predicted>0.5+pp.thres | Predicted<(0.5-pp.thres)) 
    fsam <- nrow(out)/nrow(x) # fraction of samples left 
    out <- confusionMatrix(ifelse(out$Predicted>0.5, 1, 0) %>% factor, out$Actual %>% factor)$byClass %>% data.frame
    data.frame(KIR=gsub('KIR_', '', y),
               Value=out[, 1],
               Metric=rownames(out),
               Thrs=pp.thres,
               SampleF=fsam)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)

# plot per KIR gene
p.balanced.acc.pp <- ggplot(balanced.acc.pp %>% filter(., Metric %in% c('Balanced Accuracy', 'Neg Pred Value', 'Pos Pred Value')), 
       aes(Thrs, Value, color=Metric %>% factor)) +
  geom_line() +
  xlab('PP threshold size') +
  facet_wrap(~KIR) +
  theme_minimal() +
  theme(panel.border=element_rect(fill=NA, size=0.3), legend.position='top',
        panel.grid.minor=element_blank(), legend.title=element_blank())

# fraction of samples left after PP filtering
p.sample.frac.pp <- ggplot(balanced.acc.pp %>% filter(., Metric %in% c('Balanced Accuracy', 'Neg Pred Value', 'Pos Pred Value')), 
       aes(Thrs, SampleF)) +
  geom_line() +
  xlab('PP threshold size') +
  ylab('Fraction of samples left') +
  facet_wrap(~KIR) +
  theme_minimal() +
  theme(panel.border=element_rect(fill=NA, size=0.3), legend.position='top',
        panel.grid.minor=element_blank(), legend.title=element_blank(),
        plot.margin=unit(c(44.5, 5.5, 5.5, 5.5), "pt"))



## ----------------------------------------------------------
## figures
## ----------------------------------------------------------

# Figure 2
jpeg('./results/Fig2_KIR_allplots.jpeg', width=10, height=9, res=600, units='in')
ggarrange(
  ggarrange(p.all.metrics, 
            cowplot::plot_grid(p.wgs.comp, p.kir.model.preds.missing.stats, nrow=2, ncol=1, 
                      rel_heights=c(0.45, 0.55), labels=c('b', 'c'), align='v', axis=c("l", "r")), 
            align="h",  ncol=2, labels=c('a')),#,  rel_widths=c(1, 1), axis=c("b", "t")),
  ggarrange(p.test.bee, p.test.conf, nrow=1, ncol=2, labels=c('d', 'e')),
  ncol=1, nrow=2, heights=c(0.9, 1), align='hv'
)
dev.off()


# Figure 3
jpeg('./results/Fig3_KIR_samplefrac.jpeg', width=10, height=5, res=600, units='in')
ggarrange(p.balanced.acc.pp, p.sample.frac.pp, ncol=2, labels=c('a', 'b'))
dev.off()



## ----------------------------------------------------------
## tables
## ----------------------------------------------------------

## KIR gene imputation metrics

# format metrics into columns
table.1 <- filter(all.metrics, Group=='test set') %>% dplyr::select(., -Group) %>% 
  pivot_wider(., names_from=c('Metric'), values_from='Value') %>% data.frame() %>% .[, c(1:7,12:13)]

# Table 1
fwrite(data.frame(KIR_gene=table.1[, 1], apply(table.1[, -1], 2, signif, digits=3)), './results/Table1.tsv', sep='\t')


## SNPs used in imputation

# SNPs per gene
kir.model.full.gene.var <- map2(kir.model.full, names(kir.model.full), function(x, y) {
  tmp  <- names(x$Features)
  tmp  <- str_split_fixed(tmp, '_', 3)
  data.frame(Imputation_target_gene=y, ModelVariantName=names(x$Features), Chr=19, Pos=as.numeric(tmp[, 2]), CountedAllele=tmp[, 3], 
             Importance=x$Features) %>% arrange(., 1/Importance)
}) %>% do.call(rbind, .)

# Table S1 
fwrite(kir.model.full.gene.var[, c(1, 3:6)], './results/TableS1.tsv', sep='\t')


# # get ensembl data on SNPs
# library(biomaRt)
# #ensembl <- useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")
# ensembl <- useEnsembl(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", mirror = "uswest")
# listFilters(ensembl)
# listAttributes(ensembl)
# 
# kir.model.full.vars.ensembl <- map(kir.model.full.vars$Pos, function(x) {
#   print(x)
#   getBM(c("refsnp_id", "allele", "chr_name", "chrom_start", "chrom_end", "chrom_strand"),
#         filters=c("chr_name", "start", "end"), values=list("19", x, x), mart=ensembl)
# }) %>% do.call(rbind, .)
# 
# # join with model SNPs
# table.2 <- map(kir.model.full.gene.var, function(x) {
#   left_join(x, kir.model.full.vars.ensembl, by=c('Pos'='chrom_start')) %>% 
#     dplyr::select(., -c(allele, chr_name, chrom_end, chrom_strand)) %>% .[, c(1, 7, 2:6)]
# }) %>% do.call(rbind, .)
# 
# # write out 
# fwrite(table.2, './results/Table2.tsv', sep='\t')


