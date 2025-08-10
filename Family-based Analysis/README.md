Family based analysis will be performed simply by Python and R scripts, since the GATK and slivar always fail. 

### Input files
1. Family merged vcf.gz file
2. Family ped file

### Command for analysis
```bash
newbornrisk='/mnt/storage_pool/Genomics/Genome/NewbornRisk'
bash NewbornRisk_family.sh \
  -i F4 \
  -o $newbornrisk/examples/F4 \
  -v $newbornrisk/examples/F4/F4_merged.vcf.gz \
  --ped $newbornrisk/examples/F4/F4.ped \
  --only-pass yes --genome GRCH38 --genedb TR
```
* feels like the --fork and --threads do not work, it stops the VEP processes

### UI design
#### Summary 
- Family Information: `F4_family_information.txt`
- Screen-positive Neonatal Diseases and Disease carrier status (Caused by ClinVar and ACMG (likely) pathogenic variants): `F4_plp_brief.txt`
- Variant Pedigree Analysis: pedigree figures for each variant presented in `F4_BRCA1_pedigree_plot.png`, `F4_C8B_pedigree_plot.png`, and other files ends with `pedigree_plot.png`

* here change the filename into Gene_Variant_plot.png

<img width="1920" height="968" alt="image" src="https://github.com/user-attachments/assets/e5e2181f-5fee-4816-a1d9-479c980d8f26" />
<img width="1914" height="737" alt="image" src="https://github.com/user-attachments/assets/e662948e-4134-4d15-8f65-dc333222be90" />

#### Overall Variant Statistics
- VEP statistics and figures: `F4_vep_annotated.vcf.gz_summary.html`
- Newborn Variant Statistics: `Kid0Pattern_P001_263_Variant_statistics.png`
- Loss-of-function variants and other predicted variants 
  - Loss-of-function: `F4_high_brief.txt`
  - Other predicted variants: `F4_predict_brief.txt`


#### Variant Details
`F4_detailed_table.txt`

<img width="1917" height="732" alt="image" src="https://github.com/user-attachments/assets/ff26aba1-6aff-462f-b903-83b59b456f69" />


#### PGx 
PGx results for each kids: files ends with `PGx.txt`, such as `F4_P001_263_PGx.txt`
