### Command for anlaysis
```bash
bash Scripts/FamilyRisk_cohort.sh \
	-i newborn103 \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort/newborn103_merged.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort \
	--sample-metadata /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort/newborn103gender.txt \
	--genome GRCH38 --only-pass yes --genedb NBScreening \
	--only-clinvar yes --clinvar Pathogenic,Likely_pathognic \
  	--fork 20 --threads 20 \
  	--run-prs yes --run-imputation yes --pgsid PGS002760 --run-pgx no
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
A) Table of PGS scores: `cohort_PGS.txt` \
B) Density plot compared with other populations: `PGS_Znorm1_DensityPlot.png` and `PGS_Znorm2_DensityPlot.png` \
C) Link out to the PGS-Catalog HTML file: `report.html`

### 6) Statistics of PGx summary. The reports for each individual could be downloaded as compressed file.

### 7) Polygenic risk score
