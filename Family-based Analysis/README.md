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
`family.py`
```bash
for i in {1..11}; do \
  /usr/local/bin/python3 ../Scripts/family/family.py F${i}/F${i}_merged_rmmissingalt_biallelic.vcf.gz F$i/F${i}.txt F$i/F${i}.ped;
done
```

#### 4.1 Pedigree plot in R 
`PedigreePlot.R`

#### 4.2 Rshinny for Family-based analysis page 
`family_Rshinny.R`

<img width="1920" height="968" alt="image" src="https://github.com/user-attachments/assets/e5e2181f-5fee-4816-a1d9-479c980d8f26" />
<img width="1914" height="737" alt="image" src="https://github.com/user-attachments/assets/e662948e-4134-4d15-8f65-dc333222be90" />
<img width="1917" height="732" alt="image" src="https://github.com/user-attachments/assets/ff26aba1-6aff-462f-b903-83b59b456f69" />

