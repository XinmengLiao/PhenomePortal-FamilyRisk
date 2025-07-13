Family based analysis will be performed simply by Python and R scripts, since the GATK and slivar always fail. 

#### Input files
1. Family merged vcf.gz file
2. Family ped file

#### 1. Remove missing ALT and split into biallelic
```bash
# remove missing ALT
bcftools view -e 'ALT = "."' ${family}.vcf.gz -Oz -o ${family}_rmmissingalt.vcf.gz --threads 20

# split into biallelic
bcftools norm -m -both -Oz -o ${family}_biallelic.vcf.gz ${family}_rmmissingalt.vcf.gz --threads 20
```

#### 2. VEP annotation

#### 3. Parse variants by python into text file 
`Scripts/family.py`
```bash
for i in {1..11}; do \
  /usr/local/bin/python3 ../Scripts/family/family.py F${i}/F${i}_merged_rmmissingalt_biallelic.vcf.gz F$i/F${i}.txt F$i/F${i}.ped;
done
```

#### 4. Pedigree plot in R 
`Scripts/pedigree.R`
