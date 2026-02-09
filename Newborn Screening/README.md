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
