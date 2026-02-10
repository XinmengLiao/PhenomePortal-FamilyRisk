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
The whole result page could be devided into these subpanels:
1) Cohort Information and Variant Summary: `Cohort_summary.txt`
2) Variant Details: `newborn103.txt` 
3) Variant Statistics: `plp_var_stat.txt`
4) Disease Ontology Statistics: `ontology_stat.txt`
5) Screen-positive monogenic rare disease (ClinVar **AND** ACMG P/lP variants): `monogenic_positive_results.txt`
6) Disease Carrier status (ClinVar **AND** ACMG P/lP variants): `carrier_status_results.txt`
7) Inhouse Allele frequency for all filtered variants: `inhouse_allele_frequency.txt`
* Note: variants on Chromosome X and Y are calculated in coincident with the gender file. 
8) PGS scores: (a) scores: `PGx_Scores/newborn103/score/cohort_PGS.txt`; (b) reports: `PGx_Scores/newborn103/score/report.html`; (c) Density plot: `PGx_Scores/newborn103/score/PGS_Znorm1_DensityPlot.png` and `PGx_Scores/newborn103/score/PGS_Znorm2_DensityPlot.png`.
9) PGx information: 103 individual reports in `PGx_Reports`.
