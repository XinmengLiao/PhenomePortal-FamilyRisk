### Overview
NewbornRisk location on server: `/mnt/storage_pool/Genomics/Genome/NewbornRisk`
1. Script for batch samples analysis: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/NewbornRisk_batch.sh`
2. Script for manage the VEP annotated VCF: ``
3. Script for manage the result output from python: `Batch20251203.R`
4. Example output files: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/newborn103`
5. Input file for example:
  - Merged VCF file: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/newborn103/newborn103_merged.vcf.gz`
  - SampleIDs and Genders: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/newborn103/newborn103gender.txt`

### Command for anlaysis
```bash
newbornrisk='/mnt/storage_pool/Genomics/Genome/NewbornRisk'
bash NewbornRisk_batch.sh \
  -i newborn103 \
  -v $newbornrisk/examples/P0064/newborn103_merged.vcf.gz \
  -o $newbornrisk/examples/P0064 \
  --sample-metadata $newbornrisk/examples/P0064/newborn103gender.txt \
  --genome GRCH38 --only-pass yes --fork 48 --threads 20 --genedb TR
```

## UI Design
The whole result page could be devided into 6 parts:
1) Cohort statistics and Variant statistics.
2) Result Summary: Screen-positive and Carrier status.
3) Inhouse Allele frequency for all filtered variants.
4) Variant Details (Table)
5) PGS density plot and the scores for all individual. (link out to the HTML file). 
6) Statistics of PGx summary. The reports for each individual could be downloaded as compressed file.


### 1) Cohort statistics and Variant statistics
Just show the table of: `Cohort_summary.txt`

### 2) Result Summary: Screen-positive and Carrier status
* Note: No matter what is the screening criteria, all the screen-positive results shown in the webpage refers to the diseases that caused by ClinVar AND ACMG P/lP variants. This is the default display and the barplot content of `plp_var_stat.txt` can not be changed by the users. \

The tables: \
  A) screen-positive monogenic result: `monogenic_positive_results.txt` \
  B) Carrier status results: `carrier_status_results.txt` \
  C) Table and barplot showing the distribution of disease ontology: `ontology_stat.txt`

### 3) Inhouse Allele frequency for all filtered variants
A table showing all the filtered variants, with allele frequency: `inhouse_allele_frequency.txt`
* Note: variants on Chromosome X and Y are calculated in coincident with the gender file. 

### 4) Variant Details
A table showing all the variant annotations: `newborn103_vep_merged_rmmissingalt_biallelic.txt`

### 5) PGS density plot and the scores for all individual
A) Table of PGS scores: `newborn103_pgs.txt` \
B) Density plot compared with other populations: `PRS_density.png` \
C) Link out to the PGS-Catalog HTML file: `report.html`

### 6) Statistics of PGx summary. The reports for each individual could be downloaded as compressed file.

<img width="887" height="662" alt="image" src="https://github.com/user-attachments/assets/2c52380a-966a-4a4a-85bc-02607de819f5" />
<img width="877" height="602" alt="image" src="https://github.com/user-attachments/assets/58499410-8405-4685-876b-8b3198f31d80" />
<img width="1670" height="493" alt="image" src="https://github.com/user-attachments/assets/8356c5dc-4076-4bdb-bc73-b4a254f78cea" />

