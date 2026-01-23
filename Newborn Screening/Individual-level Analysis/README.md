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

# single - newborn with PGS
# If run PGS, must avoid using underscore in naming files and folders
bash $Scripts/FamilyRisk_single.sh --newborn \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid/RapidWGS_001_C.hard-filtered.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid \
	--genome GRCH38 --only-pass yes --gender Female --genedb NBScreening --run-pgx yes \
  	--fork 20 --threads 20 \
  	--run-prs yes --run-imputation yes --pgsid PGS002760
```

## UI design
Results are stored in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid/Results`
The UI could be similar to the original newborn screening in XOmics. Functions such as selecting the variants, genes, and add comments and suggestions could be included in the web page. 
It would be better if the following details could be presented: 
1) The complete result table of filtered variants: `rwgs1_kid.txt`
2) Sample and variant summary: `general_summary.txt`
3) Screen-positive results (only caused by ClinVar **and** ACMG P/LP variants): `carrier_status_results.txt`
4) Disease carrier status (only caused by ClinVar **and** ACMG P/LP variants): `monogenic_positive_results.txt`
5) PGx information: `PGx_Reports`
6) PGS information: `PGS_Score`
