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

#### Summary 
<img width="887" height="662" alt="image" src="https://github.com/user-attachments/assets/2c52380a-966a-4a4a-85bc-02607de819f5" />
<img width="867" height="599" alt="image" src="https://github.com/user-attachments/assets/de068e9e-7def-4e05-97be-7ac8d24598f4" />
<img width="1670" height="493" alt="image" src="https://github.com/user-attachments/assets/8356c5dc-4076-4bdb-bc73-b4a254f78cea" />


#### Vairnat Details
<img width="1671" height="908" alt="image" src="https://github.com/user-attachments/assets/a88217e8-b8e3-48b3-97ec-0425e16a2371" />


