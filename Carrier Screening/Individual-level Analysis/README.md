## Carrier-Individual
### Command for analysis 
```bash
# single - carrier (without PGS) 
Scripts="/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts"
bash $Scripts/FamilyRisk_single.sh --carrier \
	-i rwgs1_mother \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-carrier/RapidWGS_001_M.hard-filtered.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-carrier \
	--genome GRCH38 --only-pass yes --gender Female --run-pgx yes \
  	--fork 20 --threads 20 
```

### UI design
Results are stored in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-carrier/Results`
The UI design could be similar to Individual-level Newborn Screening. 
1) Variant Details: `rwgs1_mother.txt`
2) Potential recessive disease that will be inherited to offsprings (ClinVar **AND** ACMG pathogenic variants): `carrier_screening_result_single.txt`
3) PGx information: `PGx_Reports/rwgs1_mother_biallelic_nodup_pass.report.html`
