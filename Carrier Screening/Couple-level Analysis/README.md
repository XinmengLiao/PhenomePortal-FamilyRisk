### Command for analysis
```bash
# family - carrier without PGS
# the gene-disease db is customized with all ACMG-carrier-screening-list and IDH3A added. 
bash Scripts/FamilyRisk_family.sh --carrier \
	-i rwgsF1 \
	--sample-list /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1/file_list.txt \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1 \
	--ped /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1/rwgs_F1.ped \
	--customized-genedb /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-carrier-rwgs1/carrier_rwgs1_customized_genelist.txt \
	--genome GRCH38 --only-pass yes --carrier --run-pgx yes \
  	--fork 20 --threads 20
```

### UI design
It will be nice if the following contents could be included in the webpage:
1) Couple details and variant summary.
2) Potential recessive disease that will be inherited to the offsprings.
3) PGS information.
4) PGx information.

#### 1) Couple details and variant summary.
1. The table of family details: `family_information.txt`
2. The table of variant summary: `variant_summary.txt`

#### 2) Potential recessive disease that will be inherited to the offspring
1. Table of potential recessive diseases: `carrier_screening_result_couple.txt`
* Note: Only diseases caused by ClinVar **AND** ACMG P/LP variants will be shown as default. 

#### 3) PGx information (have not done yet)
Results from PharmCat in HTML and/or JSON file.
