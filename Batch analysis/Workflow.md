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


#### Calculate the AF (sex chromosome excluded)
```bash
# since AF calculation might be incorrect on sex chromosomes, only conduct AF calculation on normal chromosomes
vcftools --gzvcf new_merge_ref.vcf.gz \
  --chr chr1 --chr chr2 --chr chr3 --chr chr4 --chr chr5 \
  --chr chr6 --chr chr7 --chr chr8 --chr chr9 --chr chr10 \
  --chr chr11 --chr chr12 --chr chr13 --chr chr14 --chr chr15 \
  --chr chr16 --chr chr17 --chr chr18 --chr chr19 --chr chr20 \
  --chr chr21 --chr chr22 \
  --freq --out AF.autosomes
```

#### Python and configuration file to decipher the merged files. 

#### Analysed results visualizations 
