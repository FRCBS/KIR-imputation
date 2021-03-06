## KIR-imputation
KIR gene content imputation models for 

     KIR2DL1
     KIR2DL2
     KIR2DL3
     KIR2DL5
     KIR2DP1
     KIR2DS1
     KIR2DS2
     KIR2DS3
     KIR2DS4
     KIR2DS5
     KIR3DL1
     KIR3DS1


**Published article**: Ritari J, Hyvärinen K, Partanen J and Koskela S. KIR gene content imputation from single-nucleotide polymorphisms in the Finnish population.
_PeerJ_, Jan 7, 2022 (https://peerj.com/articles/12692/)

### Dependencies and requirements
[plink 1.9](https://www.cog-genomics.org/plink/) (tested on v1.90b6.6 64-bit)

[plink 2.0](https://www.cog-genomics.org/plink/2.0/) (tested on v2.00a2.3LM 64-bit Intel)

[R](https://www.r-project.org/) (tested on v4.1.1)

[ranger](https://cran.r-project.org/web/packages/ranger/index.html) (tested on v0.13.1)

[tidyverse](https://cran.r-project.org/web/packages/-tidyverse/index.html) (tested on v1.3.1)

[data.table](https://cran.r-project.org/web/packages/data.table/index.html) (tested on v1.14.2)   

Please cite the above software and our manuscript if you use this method. 

Tested on Ubuntu 20.04.3 LTS, Intel(R) Core(TM) i7-5600U, 16GB

### Code (./src)
`functions.R` Helper functions for model fitting and evaluation

`KIR_data.R` Data preparation

`KIR_models.R` RF model fitting and testing

`KIR_plotting.R` Manuscript figures

`KIR-IMP.R` Data preparation and output handling for KIR\*IMP

`run_KIR_imputation.R` Script to run KIR imputation for plink formatted input data files

`train_models.R` Script for training models on user's own phenotype and genotype files

### Models (./models)
Fitted imputation models for the 12 KIR genes as RDS files, one file per gene. 
Contains also the associated data files needed for applying the models: _plink_allele_ref_ and _SNP_data.rds_.

### Testing (./test)
Intended for testing the functionality of the scripts and providing a reference for data input formats. Contains an artificial phenotype data file for building models on the 1000 Genomes data. 
To try out the ready-made imputation models and to test fitting new models on the 1000 Genomes data, download plink formatted genotype files (on reference build GRCh37/hg19):

```
mkdir ./test/1kG_data
wget -P ./test/1kG_data https://www.dropbox.com/s/y5wva3dcz0zdb7u/chr19_phase3.pgen.zst
wget -P ./test/1kG_data https://www.dropbox.com/s/b6km708bvenlm3j/chr19_phase3.pvar.zst
wget -P ./test/1kG_data https://www.dropbox.com/s/nhfhskyy50sqsf1/phase3_orig.psam
unzstd ./test/1kG_data/chr19_phase3.pgen.zst
unzstd ./test/1kG_data/chr19_phase3.pvar.zst
```

extract the KIR genomic region and convert to the standard plink (.bed/.bim/.fam) format:
```
plink2 --pgen ./test/1kG_data/chr19_phase3.pgen \
       --pvar ./test/1kG_data/chr19_phase3.pvar \
       --psam ./test/1kG_data/phase3_orig.psam \
       --rm-dup exclude-all \
       --set-all-var-ids @_# \
       --chr 19 --from-mb 54.4 --to-mb 55.1 \
       --max-alleles 2 \
       --make-bed \
       --out ./test/1kG_data/chr19_phase3_KIR
```

**Imputing KIRs with ready-made models**

The imputation script needs a plink .bim file as well as a genotype dosage file (.raw). This has to be in a correct allele orientation, so using the allele reference file from the /models folder is required when converting to dosage format:
```
plink --bfile ./test/1kG_data/chr19_phase3_KIR \
      --recode-allele ./models/plink_allele_ref \
      --recode A \
      --out ./test/1kG_data/chr19_phase3_KIR
```

The KIR imputation script has four arguments in this order: 
* path to the folder containing the models
* genotype dosages (.raw)
* SNPs (.bim) in plink format
* path to the output folder 
```
mkdir ./test/1kG_KIR_imputation
Rscript ./src/run_KIR_imputation.R \
        ./models \
        ./test/1kG_data/chr19_phase3_KIR.raw \
        ./test/1kG_data/chr19_phase3_KIR.bim \
        ./test/1kG_KIR_imputation
```
If the script runs technically OK, it should find only 13 SNPs shared between the models and the 1000 Genomes target data. This is expected since the models and the data are on different reference builds. The imputation result files are tab-delimited text tables containing a sample ID column and imputation posterior probabilities for each class of the input phenotype (e.g. gene presence/absence).

**Training new models**

The following script can be used to train new models on input data. Model training is not limited to KIRs, but can in principle be used for any SNP-phenotype relationship the user has. The input genotype data should be in the standard plink format (.bim & .raw). The input phenotype data (_1kG_KIR_testpheno.tsv_ in this example) is a tab-delimited text file contaning the subject IDs in the first column and (KIR) phenotypes in the following columns. The reference phenotypes are treated as class variables.

```
mkdir ./test/1kG_models
Rscript ./src/train_models.R \
        ./test/1kG_data/chr19_phase3_KIR.bim \
        ./test/1kG_data/chr19_phase3_KIR.raw \
        ./test/1kG_KIR_testpheno.tsv \
        ./test/1kG_models
```
The script saves the fitted imputation models to the given output dir along with a SNP data file, a plink allele reference file and OOB error estimates of the fitted models.

