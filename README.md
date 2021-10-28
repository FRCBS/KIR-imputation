## KIR-imputation
KIR gene content imputation models for 

     <sub>KIR2DL1</sub>
     <sub>KIR2DL2</sub>
     <sub>KIR2DL3</sub>
     <sub>KIR2DL5</sub>
     <sub>KIR2DP1</sub>
     <sub>KIR2DS1</sub>
     <sub>KIR2DS2</sub>
     <sub>KIR2DS3</sub>
     <sub>KIR2DS4</sub>
     <sub>KIR2DS5</sub>
     <sub>KIR3DL1</sub>
     <sub>KIR3DS1</sub>


*Manuscript*: Ritari J, Hyvärinen K, Partanen J and Koskela S. KIR gene content imputation from single-nucleotide polymorphisms in the Finnish population. 

#### dependencies and requirements
[plink 1.9](https://www.cog-genomics.org/plink/) (tested on v1.90b6.6 64-bit)

[plink 2.0](https://www.cog-genomics.org/plink/2.0/) (tested on v2.00a2.3LM 64-bit Intel)

[R](https://www.r-project.org/) (tested on v4.1.1)

[ranger](https://cran.r-project.org/web/packages/ranger/index.html) (tested on v0.13.1)

[tidyverse](https://cran.r-project.org/web/packages/-tidyverse/index.html) tested on (v1.3.1)

[data.table](https://cran.r-project.org/web/packages/data.table/index.html) tested on (v1.14.2)   

Please cite the above software if you use this method.

Tested on Ubuntu 20.04.3 LTS, Intel(R) Core(TM) i7-5600U, 16GB

#### code (./src)
`functions.R` Helper functions for model fitting and evaluation

`KIR_data.R` Data preparation

`KIR_models.R` RF model fitting and testing

`KIR_plotting.R` Manuscript figures

`KIR-IMP.R` Data preparation and output handling for KIR\*IMP

`run_KIR_imputation.R` Script to run KIR imputation for plink formatted input data files

`train_models.R` Script for training models on user's own phenotype and genotype files

#### models (./models)
Fitted imputation models for the 12 KIR genes as .rds files, one file per gene. 
Contains also the associated data files needed for applying the models: \plink_allele_ref and \SNP_data.rds.

#### tests (./models/test)
Folder for testing and demonstration purposes. Contains an artificial phenotype data file for building models on the 1000 Genomes data. 
To try out the ready-made imputation models on the 1000 Genomes data, download plink formatted genotype files (run from the repository root)

```
mkdir ./test/1kG_data
wget -P ./test/1kG_data https://www.dropbox.com/s/y5wva3dcz0zdb7u/chr19_phase3.pgen.zst
wget -P ./test/1kG_data https://www.dropbox.com/s/b6km708bvenlm3j/chr19_phase3.pvar.zst
wget -P ./test/1kG_data https://www.dropbox.com/s/nhfhskyy50sqsf1/phase3_orig.psam
unzstd ./test/1kG_data/chr19_phase3.pgen.zst
unzstd ./test/1kG_data/chr19_phase3.pvar.zst
```

and extract the KIR region into the correct plink format
```
plink2 --pgen ./test/1kG_data/chr19_phase3.pgen \
       --pvar ./test/1kG_data/chr19_phase3.pvar \
       --psam ./test/1kG_data/phase3_orig.psam \
       --rm-dup exclude-all \
       --set-all-var-ids @_# \
       --chr 19 --from-mb 54.4 --to-mb 55.0 \
       --max-alleles 2 --make-bed \
       --out ./test/1kG_data/chr19_phase3_KIR
```

Make a dosage file (.raw) using the allele orientation reference from the models folder
```
plink --bfile ./test/1kG_data/chr19_phase3_KIR \
      --recode-allele ./models/plink_allele_ref \
      --recode A \
      --out ./test/1kG_data/chr19_phase3_KIR
```

Run KIR imputation. 
```
mkdir ./test/1kG_KIR_imputation
Rscript ./src/run_KIR_imputation.R \
        ./models \
        ./test/1kG_data/chr19_phase3_KIR.raw \
        ./test/1kG_data/chr19_phase3_KIR.bim \
        ./test/1kG_KIR_imputation
```
The output result files are tab-delimited tables containing sample ID, imputation posterior probabilities for each class of the .


An example of training models on the 1000 Genomes data. The above downloaded genotype data (.bim, .raw) are in appropriate format for fitting models. The reference phenotype data \1kG_KIR_testpheno.tsv is a tab-delimited text file contaning subject IDs in the first column and (KIR) genotypes in the following columns. The reference genotypes are treated as class variables.

```
mkdir ./test/1kG_models
Rscript ./src/train_models.R \
        ./test/1kG_data/chr19_phase3_KIR.bim \
        ./test/1kG_data/chr19_phase3_KIR.raw \
        ./test/1kG_KIR_testpheno.tsv 
        ./test/1kG_models
```
The script produces the fitted imputation models in the given output dir along with associated files for SNP data, plink allele reference and OOB error estimates.

