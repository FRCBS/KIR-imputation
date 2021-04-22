# KIR-imputation
KIR gene content imputation from FinnGen genotype platform SNP data in the Finnish population.
Manuscript: "KIR gene content imputation from single-nucleotide polymorphisms in the Finnish population". 
KIR gene absence/presence level imputation for 
KIR2DL1 - KIR2DL2 - KIR2DL3 - KIR2DL5 - KIR2DP1 - KIR2DS1 - KIR2DS2 - KIR2DS3 - KIR2DS4 - KIR2DS5 - KIR3DL1 - KIR3DS1
------------ 
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

#### code (./src)

#### results (./results)
Figures and tables in the manuscript.


#### models (./models)
Fitted models for the 12 KIR genes as .rds files.

#### test data (./test)
The test data folder contains artificial genotype data that can be used as a reference when aligning SNP orientations of an input dataset.
For example, [GenotypeHarmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) by [Deelen et al.](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-7-901) can be used to orient new input data as follows:
```shell
    java -Xmx1g -jar GenotypeHarmonizer.jar \
        --input ./data/mydata \
        --inputType PLINK_BED \
        --ref ./data/test/simulated_ref \ 
        --refType PLINK_BED \
        --output ./data/mydata_harmonized \
        --outputType PLINK_BED \ 
        --update-id --keep --update-reference-allele
    
    plink --bfile ./data/mydata_harmonized --keep-allele-order --recode A --out ./data/mydata_harmonized_dosages
```   
The output allele dosage data can then be used as an input to the KIR imputation pipeline.

