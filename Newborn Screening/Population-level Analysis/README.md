## Newborn-Population
### Command for anlaysis
```bash
# for 103 TR newborns
Scripts="/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts"
bash $Scripts/FamilyRisk_cohort.sh \
	-i newborn103 \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort/newborn103_merged.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort \
	--sample-metadata /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort/newborn103gender.txt \
	--genome GRCH38 --only-pass yes --genedb NBScreening \
	--only-clinvar yes --clinvar Pathogenic,Likely_pathognic \
  	--fork 20 --threads 20 \
  	--run-prs yes --run-imputation yes --pgsid PGS002760 --run-pgx no

# smaller test cohort, four newborns were extracted from the merged file. 
bash Scripts/FamilyRisk_cohort1.sh \
	-i cohort_small \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort/cohort_small_biallelic_nodup_pass.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort \
	--sample-metadata /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort/cohort_small_gender.txt \
	--genome GRCH38 --only-pass yes --genedb NBScreening \
	--only-clinvar yes --clinvar Pathogenic,Likely_pathognic --af-clinvar 0.05 \
  	--fork 20 --threads 20 \
  	--run-prs yes --run-imputation yes --pgsid PGS002248 --run-pgx no
```


### UI Design
Results are stored in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/TRcohort/Results`
The whole result page could be devided into 6 parts:
1) Complete results of filtered variants.
2) Cohort statistics and Variant statistics.
3) Result Summary: Screen-positive and Carrier status.
4) Inhouse Allele frequency for all filtered variants.
5) PGS density plot and the scores for all individual. (link out to the HTML file).

#### 1) Complete result table.
`newborn103.txt` 

#### 2) Cohort statistics and Variant statistics
Just show the table of: `Cohort_summary.txt`

#### 3) Result Summary: Screen-positive and Carrier status
* Note: No matter what is the screening criteria, all the screen-positive results shown in the webpage refers to the diseases that caused by ClinVar AND ACMG P/lP variants. This is the default display and the barplot content of `plp_var_stat.txt` can not be changed by the users. \

The tables: \
  A) screen-positive monogenic result: `monogenic_positive_results.txt` \
  B) Carrier status results: `carrier_status_results.txt` \
  C) Table and barplot showing the distribution of disease ontology: `ontology_stat.txt`

#### 4) Inhouse Allele frequency for all filtered variants
A table showing all the filtered variants, with allele frequency: `inhouse_allele_frequency.txt`
* Note: variants on Chromosome X and Y are calculated in coincident with the gender file. 

#### 5) PGS density plot and the scores for all individual
A) Table of PGS scores: `cohort_PGS.txt` \
B) Density plot compared with other populations: `PGS_Znorm1_DensityPlot.png` and `PGS_Znorm2_DensityPlot.png` \
C) Link out to the PGS-Catalog HTML file: `report.html`
