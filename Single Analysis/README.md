## Overview 
NewbornRisk location on server: `/mnt/storage_pool/Genomics/Genome/NewbornRisk`
1. Script for single analysis: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/NewbornRisk_single.sh`
2. Script for manage the VEP annotated VCF: ``
3. Script for manage the result output from python: `Single20251203.R`
4. Example output files: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/single`
5. Input file for example:
  - Merged VCF file: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/single/P0064_1203.vcf.gz`

## Command for analysis 
```bash
newbornrisk='/mnt/storage_pool/Genomics/Genome/NewbornRisk'
bash NewbornRisk_single.sh \
  -i P0064_1203 \
  -o $newbornrisk/examples/single \
  -v $newbornrisk/examples/single/P0064_1203.vcf.gz \
  --genome GRCH38 --only-pass yes --gender Female --genedb TR \
  --fork 20 --threads 20
```

## Output files
- VEP annotated file: `P0064_1203_vep_annotated.vcf.gz`
- Detailed nodup table: `P0064_1203.txt` or `rwgs_001_C.txt`

## UI design
The UI could be similar to the original newborn screening in XOmics. Functions such as selecting the variants, genes, and add comments and suggestions could be included in the web page. 
It would be better if the following details could be presented: 
1) Sample and variant summary: `general_summary.txt`
2) Screen-positive results (only caused by ClinVar and ACMG P/LP variants): `carrier_status_results.txt`
3) Disease carrier status (only caused by ClinVar and ACMG P/LP variants): `monogenic_positive_results.txt`
4) PGx information (HTML/JSON from PharmCat)
