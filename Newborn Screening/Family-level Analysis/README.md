### Overview
NewbornRisk location on server: `/mnt/storage_pool/Genomics/Genome/NewbornRisk`
1. Script for family analysis: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/NewbornRisk_family.sh`
2. Script for manage the annotated variants from VEP: ``
3. Script for manage the python result: `Pedigree_Analysis20251203.R`
4. Example output files: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/F4`
5. Input file for example:
   - Merged VCF file: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/F4/F4_merged.vcf.gz`
   - Family information PED file: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/F4/F4.ped`
  
### Command for analysis
```bash
newbornrisk='/mnt/storage_pool/Genomics/Genome/NewbornRisk'
bash NewbornRisk_family.sh \
  -i F4 \
  -o $newbornrisk/examples/F4 \
  -v $newbornrisk/examples/F4/F4_merged.vcf.gz \
  --ped $newbornrisk/examples/F4/F4.ped \
  --only-pass yes --genome GRCH38 --genedb TR \
  --fork 20 --threads 20
```

### UI design
It will be nice if the following contents could be included in the webpage:
1) Family details and variant summary.
2) Screen-positive rare diseases and carrier status.
3) Pedigree plots for the variants of screen-positive diseases.
4) Pedigree plot generator for varant of interest.
5) PGS information.
6) PGx information.

#### 1) Family details and variant summary.
1. The table of family details: `family_information.txt`
2. The table of variant summary: `variant_summary.txt`

#### 2) Screen-positive rare diseases and carrier status
1. Table of screen-positive rare disease: `carrier_status_results.txt`
2. Table of carrier status: `monogenic_positive_results.txt`
* Note: Only diseases caused by ClinVar or ACMG P/LP variants will be shown as default. 

#### 3) Pedigree plots for the variants of screen-positive diseases.
The plot will be automatical generated: `F1_IDH3A_c.802G>A_pedigree_plot.png`

#### 4) Pedigree plot generator for varant of interest
Users could choose their variant of interest to generate additional pedigree plot showing the inheritance paths. 
Script for this function: `Pedigree_Analysis_subset20251203.R`

#### 5) PGS information
1. Table of the family PGS scores: `cohort_PGS.txt`
2. Density plot compared with the reference populations: `PGS_Znorm1_DensityPlot.png` and `PGS_Znorm2_DensityPlot.png`
3. Report for the PGS score: `report.html`

#### 6) PGx information
Results from PharmCat in HTML and/or JSON file. Have not done yet. 

<img width="1920" height="968" alt="image" src="https://github.com/user-attachments/assets/e5e2181f-5fee-4816-a1d9-479c980d8f26" />
<img width="1914" height="737" alt="image" src="https://github.com/user-attachments/assets/e662948e-4134-4d15-8f65-dc333222be90" />
<img width="1917" height="732" alt="image" src="https://github.com/user-attachments/assets/ff26aba1-6aff-462f-b903-83b59b456f69" />
