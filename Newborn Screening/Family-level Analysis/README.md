## Newborn-Family
### Command for analysis
```bash
# trio data
Scripts="/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts"
bash $Scripts/FamilyRisk_family.sh --newborn \
	-i rwgsF1 \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples/rwgsF1/rwgs_F1.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-newborn-rwgs1 \
	--ped /mnt/nas/Genomics/Genome/FamilyRisk/examples/rwgsF1/rwgs_F1.ped \
	--genome GRCH38 --only-pass yes \
	--customized-genedb /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1/carrier_rwgs1_customized_genelist.txt \
  	--fork 20 --threads 20 \
  	--run-prs yes --run-imputation yes --pgsid PGS002248 --run-pgx yes
```

### UI design
Results are stored in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-newborn-rwgs1/Results`
It will be nice if the following contents could be included in the webpage:
1) Family details: `family_information.txt`
2) Variant summary: `variant_summary.txt`
3) Variant Details: `rwgsF1.txt`
4) Screen-positive monogenic rare diseases (ClinVar **AND** ACMG P/LP variants): `monogenic_positive_results.txt`
5) Disease carrier status (ClinVar **AND** ACMG P/LP variants): `carrier_status_results.txt`
6) Pedigree plots for the variants of screen-positive diseases: Nothing for this example. 
7) Pedigree plot generator for varant of interest: `F1_IDH3A_c.802G>A_pedigree_plot.png`
Users could choose their variant of interest to generate additional pedigree plot showing the inheritance paths. 
Script for this function: `/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts/RScripts/Pedigree_Analysis_subset20260121.R`
8) PGx report: Every single html file. `PGx_Reports/rwgsF1_biallelic_nodup_pass.NW_001_C.report.html`, `PGx_Reports/rwgsF1_biallelic_nodup_pass.NW_001_F.report.html`, `PGx_Reports/rwgsF1_biallelic_nodup_pass.NW_001_M.report.html`.
9) PGS information: (a) scores: `PRS_Scores/cohort.txt`; (b) report: `PGS_Scores/report.html`; (c) Densityplots: `PRS_Scores/PGS_Znorm1_DensityPlot.png` and `PRS_Scores/PGS_Znorm2_DensityPlot.png`