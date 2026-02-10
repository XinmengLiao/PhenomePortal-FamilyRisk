## Newbornx-Individual
### Command for analysis 
```bash
Scripts="/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts"
bash $Scripts/FamilyRisk_single.sh --newborn \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid/RapidWGS_001_C.hard-filtered.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid \
	--genome GRCH38 --only-pass yes --gender Female --genedb NBScreening --run-pgx yes \
  	--fork 20 --threads 20

# single - newborn with PGS (rWGS-F1-Child)
# If run PGS, must avoid using underscore in naming files and folders
bash $Scripts/FamilyRisk_single.sh --newborn \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid/RapidWGS_001_C.hard-filtered.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid \
	--genome GRCH38 --only-pass yes --gender Female --genedb NBScreening --run-pgx yes \
  	--fork 20 --threads 20 \
  	--run-prs yes --run-imputation yes --pgsid PGS002760

# single - newborn with PGS (P0064_1002)
bash Scripts/FamilyRisk_single.sh --newborn \
	-i P0064_1002 \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-1002/P0064_1002.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-1002 \
	--genome GRCH38 --only-pass yes --gender Male --genedb NBScreening --af-clinvar 0.05 --only-clinvar yes --clinvar Pathogenic,Likely_pathognic \
  	--fork 20 --threads 20 \
  	--run-prs yes --run-imputation yes --pgsid PGS001893 --run-pgx yes
```

## UI design
#### P0064_1002 is the example. 
Results are stored in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-1002/Results`
The UI could be similar to VarXOmics. Each table could be in a separated subpanel. The reporting system could be simiarl to XOmics report generating function. Newborn information, positive screening results, customized screening results, and comments & suggestions will be included in the report.  
Tables for subpanels：
1) Sample and variant summary: `general_summary.txt`
2) Variant Details: `P0064_1002.txt`
3) Disease carrier status (only caused by ClinVar **and** ACMG P/LP variants): `monogenic_positive_results.txt`
4) Screen-positive results (only caused by ClinVar **and** ACMG P/LP variants): `carrier_status_results.txt`
5) PGx information: `PGx_Reports/P0064_1002_biallelic_nodup_pass.report.html`
6) PGS information: (a) scores: `P0064-1002_pgs.txt` (need to decompressed the original file) ; (b) report: `PGS_Score/score/report.html`
