## Newborn-Family
### Command for analysis
```bash
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
1) Complete results of filtered variants. 
2) Family details and variant summary.
3) Screen-positive rare diseases and carrier status.
4) Pedigree plots for the variants of screen-positive diseases. 
5) Pedigree plot generator for varant of interest.
6) PGx information: `PGS_Reports/`
7) PGS information.

#### 1) Complete result table.
`rwgsF1.txt`

#### 2) Family details and variant summary.
1. The table of family details: `family_information.txt`
2. The table of variant summary: `variant_summary.txt`

#### 3) Screen-positive rare diseases and carrier status
1. Table of screen-positive rare disease: `carrier_status_results.txt`
2. Table of carrier status: `monogenic_positive_results.txt`
* Note: Only diseases caused by ClinVar **AND** ACMG P/LP variants will be shown as default. 

#### 4) Pedigree plots for the variants of screen-positive diseases.
The plot will be automatical generated: `F1_IDH3A_c.802G>A_pedigree_plot.png`

#### 5) Pedigree plot generator for varant of interest
Users could choose their variant of interest to generate additional pedigree plot showing the inheritance paths. 
Script for this function:

#### 6) PGx information
Results from PharmCat in HTML and/or JSON file. Have not done yet. 

#### 7) PGS information
1. Table of the family PGS scores: `cohort_PGS.txt`
2. Density plot compared with the reference populations: `PGS_Znorm1_DensityPlot.png` and `PGS_Znorm2_DensityPlot.png`
3. Report for the PGS score: `report.html`