## Carrier-Family
### Command for analysis
```bash
# family - carrier with PGS
# the gene-disease db is customized with all ACMG-carrier-screening-list and IDH3A added. 
Scripts="/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts"
bash $Scripts/FamilyRisk_family.sh --carrier \
	-i rwgsF1 \
	--sample-list /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1/file_list.txt \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1 \
	--ped /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1/rwgs_F1.ped \
	--customized-genedb /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1/carrier_rwgs1_customized_genelist.txt \
	--genome GRCH38 --only-pass yes --carrier --run-pgx yes --af-clinvar 0.05 \
	--fork 20 --threads 20 --pgsid PGS002106
```

### UI design
Results are stored in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1/Results`
It will be nice if the following contents could be included in the webpage:
1) Couple Details: `family_information.txt`
2) Variant Summary: `variant_summary.txt`
3) Variant Details: `rwgsF1.txt`
4) Potential recessive disease that will be inherited to the offspring (ClinVar **AND** ACMG P/LP variants):  `carrier_screening_result_couple.txt`
4) PGx information: report for each individual in `PGx_Reports/`
5) PGS information: (a) scores: `PGS_Scores/cohort_PGS.txt`; (b) report: `PGS_Scores/report.html`; (c) Density plots: `PGS_Znorm1_DensityPlot.png` and `PGS_Znorm2_DensityPlot.png`
