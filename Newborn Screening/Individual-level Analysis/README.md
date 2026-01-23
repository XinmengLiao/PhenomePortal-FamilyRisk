## Newborn - Individual

## Command for analysis 
```bash
# single - newborn 
bash Scripts/FamilyRisk_single.sh --newborn \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid/RapidWGS_001_C.hard-filtered.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid \
	--genome GRCH38 --only-pass yes --gender Female --genedb NBScreening --run-pgx yes \
  	--fork 20 --threads 20

# single - newborn with PGS
# If run PGS, must avoid using underscore in naming files and folders
bash Scripts/FamilyRisk_single.sh --newborn \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid/RapidWGS_001_C.hard-filtered.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-newborn-rwgs1Kid \
	--genome GRCH38 --only-pass yes --gender Female --genedb NBScreening --run-pgx yes \
  	--fork 20 --threads 20 \
  	--run-prs yes --run-imputation yes --pgsid PGS002760
```

## Output files
- VEP annotated file: `P0064_1203_vep_annotated.vcf.gz`
- Detailed nodup table: `rwgs_001_C.txt`

## UI design
The UI could be similar to the original newborn screening in XOmics. Functions such as selecting the variants, genes, and add comments and suggestions could be included in the web page. 
It would be better if the following details could be presented: 
1) Sample and variant summary: `general_summary.txt`
2) Screen-positive results (only caused by ClinVar **and** ACMG P/LP variants): `carrier_status_results.txt`
3) Disease carrier status (only caused by ClinVar **and** ACMG P/LP variants): `monogenic_positive_results.txt`
4) PGx information (HTML/JSON from PharmCat)
