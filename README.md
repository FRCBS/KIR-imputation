# KIR-imputation
KIR gene content imputation from FinnGen genotype platform SNP data in the Finnish population

### code (./src)

### test data (./test)

The test data folder contains artificial genotype data that can be used as a reference when aligning SNP orientations of an input dataset.
For example, [GenotypeHarmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) by [Deelen et al.](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-7-901) can be used to orient new input data as follows:

    java -Xmx1g -jar GenotypeHarmonizer.jar \
    --input ./data/mydata --inputType PLINK_BED \
    --ref ./data/test/simulated_ref --refType PLINK_BED \
    --output ./data/mydata_harmonized --outputType PLINK_BED \ 
    --update-id --keep --update-reference-allele
    
    plink --bfile ./data/mydata_harmonized --keep-allele-order --recode A --out ./data/mydata_harmonized_dosages
   
The output allele dosage data can then be used as an input to the KIR imputation pipeline.

