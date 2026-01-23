## Command for analysis 
```bash
# single - carrier (without PGS) 
bash Scripts/FamilyRisk_single.sh --carrier \
	-i rwgs1_mother \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-carrier/RapidWGS_001_M.hard-filtered.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single-carrier \
	--genome GRCH38 --only-pass yes --gender Female --run-pgx yes \
  	--fork 20 --threads 20 
```

## Output files
- VEP annotated file: `P0064_1203_vep_annotated.vcf.gz`
- Detailed nodup table: `P0064_1203.txt` or `rwgs_001_C.txt`

## UI design
The UI could be similar to the original newborn screening in XOmics. Functions such as selecting the variants, genes, and add comments and suggestions could be included in the web page. 
It would be better if the following details could be presented: 
1) Potential recessive disease that will be inherited to offsprings: `carrier_screening_result_single.txt`
2) PGx information 
