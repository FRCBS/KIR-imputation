## KIR-imputation
KIR gene content imputation models for 
* <sub>KIR2DL1</sub>
* <sub>KIR2DL2</sub>
* <sub>KIR2DL3</sub>
* <sub>KIR2DL5</sub>
* <sub>KIR2DP1</sub>
* <sub>KIR2DS1</sub>
* <sub>KIR2DS2</sub>
* <sub>KIR2DS3</sub>
* <sub>KIR2DS4</sub>
* <sub>KIR2DS5</sub>
* <sub>KIR3DL1</sub>
* <sub>KIR3DS1</sub>

Manuscript: Ritari J, Hyvärinen K, Partanen J and Koskela S. KIR gene content imputation from single-nucleotide polymorphisms in the Finnish population. 


#### code (./src)
`functions.R` Helper functions for model fitting and evaluation

`KIR_data.R` Data preparation

`KIR_models.R` RF model fitting and testing

`KIR_plotting.R` Manuscript figures

`KIR-IMP.R` Data preparation and output handling for KIR\*IMP

`run_KIR_imputation.R` Script to run KIR imputation for plink formatted input data files


#### models (./models)
Fitted models for the 12 KIR genes as .rds files.

##### test data (./models/test)
The test data folder contains artificial genotype data for checking the scripts. The `plink_allele_ref` file is an allele orientation reference for conversion to dosage format, for example:
```
plink --bfile ./test/simulated_ref \
      --recode-allele ./test/plink_allele_ref \
      --recode A \
      --out ./test/simulated_ref
```   
The variant names in the input data should be formatted as chr19_hg38position_allele1_allele2 , for example `chr19_54480203_A_G`.
Correctly oriented allele dosage data can then be used as an input to the KIR imputation pipeline. An example using the test files:
```
Rscript ./src/run_KIR_imputation.R \
    "./results/models" \
    "./test/simulated_ref.raw" \
    "./test/simulated_ref.bim" \
    "./test/output"
```
The output data files are tab delimited tables containing sample ID, imputation result (1 or 0 for presence/absence), and posterior probability for gene presence.


