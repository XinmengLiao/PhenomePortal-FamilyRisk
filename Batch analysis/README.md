Here is the workflow commands for batch analysis of the 103 Turkish newborns recruited from APMI. 

#### Merged batch sample preparation
##### Method 1: merge vcf files and VEP annote next 
```bash
# split multiallelic into biallelic
bcftools norm -m -both -Oz -o ${sampleID}.biallelic.vcf.gz ${sample}.vcf.gz --threads 10

# index
bcftools index ${sampleID}.biallelic.vcf.gz --threads 10

# merge (here the missing alleles are still ./.)
bcftools merge -l nbrisk.txt -Oz -o newborn103_merged.vcf.gz --threads 10 -W

# change missing allele into reference allele
bcftools +setGT newborn103_merged.vcf.gz -- -t . -n 0/0 | bgzip > newborn103_merged_ref.vcf.gz

# Annotated batch sample with VEP v113
```

##### Method 2: annotate VEP first and merge next 
```bash
# split multiallelic into biallelic
bcftools norm -m -both -Oz -o ${sampleID}.vep.biallelic.vcf.gz ${sample}.vep_annotated.vcf.gz --threads 10

# index
bcftools index ${sampleID}.vep.biallelic.vcf.gz --threads 10

# merge (here the missing alleles are still ./.)
bcftools merge -l nbrisk.txt -Oz -o newborn103_vep_merged.vcf.gz --threads 10 -W

# change missing allele into reference allele
bcftools +setGT newborn103_vep_merged.vcf.gz -- -t . -n 0/0 | bgzip > newborn103_vep_merged_ref.vcf.gz
```

#### Quality control (remove missing ALT)
`bcftools view -e 'ALT = "."' newborn103_vep_merged.vcf.gz -Oz -o newborn103_vep_merged_rmmissingalt.vcf.gz --threads 10`


#### Matching witht the gene-disease list 


#### Python and configuration file to decipher the merged files. 

#### Analysed results visualizations 
