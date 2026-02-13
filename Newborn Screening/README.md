Log: 20260213

### Updates needed for more functions
1. For all modules that generate PDF report, add a function where users could select any variant of interest from "Variant Details" and added it to the final PDF report.
2. Add the function to generate pedigree plot for any variant of interest from "Variant Details" in the family-level carrier/newborn screening analysis.
   - R script for generating additional pedigree plot: `/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts/RScripts/Pedigree_Analysis_subset20260121.R`. Here, the input is txt/tsv file of selected variant separted by line. Input file example could be found in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-newborn-rwgs1/variant_info.txt`. Every variant will have one plot generated in the output directory.
   - Command for running this Rscript: `Rscript Pedigree_Analysis_subset20260121.R $INPUT_SAMPLE $OUTPUT_DIR/${INPUT_SAMPLE}.txt $PED $OUTPUT_DIR $VAR_FILE`. Please also see the arguments in Rscript. These variables are same as `FamilyRisk_family.sh`, only $VAR_FILE is newly defined.
   - The story for this function is: In the family-level newborn screening on the trio family, we applied the customized gene-disease screening list where IDH3A-Mitochondrial encephalopathy is included (`/mnt/nas/Genomics/Genome/FamilyRisk/examples/rwgsF1/customized_genedb.txt`). In the output results, though none of the rare monogenic diseases are screened positively, clinicians find out the IDH3A c.802G>A p.Gly268Arg (chr15_78165014_G_A) may be the pathogenic variant leading to clinical phenotypes. Thus, clinician manually selected this variant from "Variant Details" and generated extra pedigree plot for this variant. Eventually, this variant is selected to be included in the final pdf report together with the pedigree plot. None of the disease carrier status was found.
   - Plot for this IDH3A variant is in: `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/family-newborn-rwgs1/Results/F1_IDH3A_c.802G>A_pedigree_plot.png`
3. Move FamilyRisk Database from Scilife server to our own server. Src in docker: `xinmengliao/familyriskdatabase` or `docker run xinmengliao/familyriskdatabase:v1`. 
4. Allow email for login/register to store personal data?

### Extra work
1. Use FamilyRisk to conduct carrier and newborn screening for two more families to test the pipeline and see if it could provide true results. 

<br>
<br>

Log: 20260211

### Updates needed for results interface
1.For PGx and PGS html reports, I can not jump to any hyperlink within the report. Error shows "You need to enable JavaScript to run this app." \
<img width="383" height="89" alt="image" src="https://github.com/user-attachments/assets/0a3924ff-7f89-49a9-89df-992ceaf08e7f" /> \
<br>
2. For PGS section, we don't need to show the json and version files. Keep the other files. \
<img width="350" height="331" alt="image" src="https://github.com/user-attachments/assets/aacc0fb2-d897-452e-a7f6-40b622272f76" /> \
<br>
3. We don't need to show "Other Files". \
<img width="344" height="126" alt="image" src="https://github.com/user-attachments/assets/ba429e21-7a0f-40e4-ba91-dfa9e538b15c" />


<br>
<br>


Log: 20260209

### Roadmap for endpoint pipeline
<img width="2456" height="1961" alt="Picture 1" src="https://github.com/user-attachments/assets/407f22b5-f73d-451b-adbd-d6cc5f200221" />

### Updates needed on the interfaces
1. Updates Gene-disease panel list with the following names:
   - Newborn screening panels: ACMG Secondary Findings (ACMGv3.3), BabySeq_GroupA, BabySeq_GroupB, EarlyCheck_Group2, EarlyCheck_Group1, EarlyCheck_Group3, EarlyCheck_Group4, Genomic101, NBScreening, BabyScreen+, BabyDetect, Guardian_Group1, Guardian_Group2
   - ACMG Expanded Carrier Screening list (ACMG-ECS): ACMG_Carrier_Tier_1, ACMG_Carrier_Tier_2, ACMG_Carrier_Tier_3, ACMG_Carrier_Tier_4, ACMG_Carrier_Outside_gnomAD
   - * Note: parameter for `--genedb` for the `ACMG Secondary Findings (ACMGv3.3)` is `--genedb ACMGv3.3`. The other parameters are same as the name listing here. 
2. Could we move CNV/SV under Newborn Screening section. Since this function is currently developed for Newborn Screening only.
3. ClinVar classifications has **Benign** and **Likely benign** as well. Please include these two in as options. But these two are not the default setting.
4. In the description for "VCF file" of Individual-level Analysis, please change "bgzipped recommended" to "bgzipped required".
5. In the description for "Gene-disease panel" of Individual-level Analysis, please change " If omitted, the script default will be used." to "If omitted, NBScreening panel will be used for newborn screening, and ACMG-ECS will be used for carrier screening."
6. In the description for "Run genotype imputation" of Individual-level Analysis, please change "Recommended for WGS with 0/0 genotypes." to "Recommended for WGS with missing genotypes."
7. In the description for "Run genotype imputation" of Family-level Newborn screening, instead of the second VCF file, users should be able to upload more than two VCF files for merging, since some families have more than three members. Could we change the function to 1) upload merged family VCF; or 2) upload multiple unmerged VCF files.
8. If we do not need JSON file for PGx, we can remove `-reporterJson` in pharmcat step B in each bash script. 
9. We need to design the reports for individual- and family-level carrier and newborn screenings.
10. Add one more function for Family-level newborn screening, where users could select specific variants and create the pedigree plots. Scripts are available at ``.

### Interface and Reports for Demos 
1. Individual-level newborn screening (SNVs and INDELs) + Reporting: https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/tree/main/Newborn%20Screening/Individual-level%20Analysis
2. Individual-level newborn screening (CNVs and SVs):https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/tree/main/Newborn%20Screening/Struactural-Variant%20Analysis
3. Family-level newborn sceening + Reporting:https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/tree/main/Newborn%20Screening/Family-level%20Analysis
4. Cohort-level newborn screening: https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/tree/main/Newborn%20Screening/Population-level%20Analysis
5. Family-level carrier screening: https://github.com/XinmengLiao/PhenomePortal-FamilyRisk/tree/main/Carrier%20Screening/Couple-level%20Analysis

### Reporting (need to design later)
1. Individual-level newborn screening
2. Individual-level carrier screening 
3. Family-level newborn screening 
4. Family-level carrier screening 
