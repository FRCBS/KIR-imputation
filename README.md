## KIR-imputation
KIR gene content imputation models for: 
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

#### results (./results)
Manuscript figures and tables.


#### models (./models)
Fitted models for the 12 KIR genes as .rds files.

#### test data (./test)
The test data folder contains artificial genotype data that can be used as a reference when aligning SNP orientations of an input dataset.
For example, [GenotypeHarmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) by [Deelen et al.](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-7-901) can be used to orient new input data as follows:
```
java -Xmx1g -jar GenotypeHarmonizer.jar \
    --input ./data/mydata \
    --inputType PLINK_BED \
    --ref ./data/test/simulated_ref \ 
    --refType PLINK_BED \
    --output ./data/mydata_harmonized \
    --outputType PLINK_BED \ 
    --update-id --keep --update-reference-allele
 
plink --bfile ./data/mydata_harmonized \
    --keep-allele-order --recode A \
    --out ./data/mydata_harmonized_dosages
```   
The output allele dosage data can then be used as an input to the KIR imputation pipeline.

