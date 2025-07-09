Here is the workflow commands for batch analysis of the 103 Turkish newborns recruited from APMI. 

#### Merged batch sample preparation
```bash
# split multiallelic into biallelic
bcftools norm -m -both -Oz -o ${sampleID}.biallelic.vcf.gz ${sample}.vcf.gz --threads 10

# index
bcftools index ${sampleID}.biallelic.vcf.gz --threads 10

# merge (here the missing alleles are still ./.)
bcftools merge -l nbrisk.txt -Oz -o newborn103_merged.vcf.gz --threads 10 -W

# change missing allele into reference allele
bcftools +setGT newborn103_merged.vcf.gz -- -t . -n 0/0 -Oz -o newborn103_merged_ref.vcf.gz

```

#### Annotated batch sample with VEP v113
```bash
# need to change ClinVar into the complete ones
```

#### Python and configuration file to decipher the merged files. 

#### Analysed results visualizations 
