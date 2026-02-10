### Commands 
```bash
# single - cnv/sv with snv
Scripts="/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts"
bash $Scripts/FamilyRisk_single_cnvsv.sh \
	-i P0064_1002 \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv/P0064_1002.sv.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv \
	--genome GRCH38 \
	--snvindel-vcf /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_newborn/RapidWGS_001_C.hard-filtered.vcf.gz

# single - cnv/sv without snv
bash $Scripts/FamilyRisk_single_cnvsv.sh \
	-i P0064_1002 \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv/P0064_1002.cnv.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv \
	--genome GRCH38 --data-type cnv
```

### UI design
Results are stored in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv/Results`
There is two types of results: 1) with compound heterozygosity; 2) without compound heterozygosity. \
If users provide `--snvindel-vcf`, then the results will be with compound heterozygositi, where the filename ends with `.cnvsv.with_compound_het.tsv`. Otherwise, the filename ends with `annotated.tsv`.  \
All CNV/SV will be filtered by only keeping the genes that are included in the screening list. So users are allowed to set `--genedb` and `--customized-genedb`. If both parameters are empty, the default NBScreening list will be applied. (all these have been already set in Python script).  \
All results will be filtered with matching genes included in the screening panel. 

1) Complete results of annotated CNV results: `RapidWGS_001_C.cnv.annotated.tsv`
2) Complete reullts of annotated CNV/SV with compound heterozygous variants identified with SNV/INDEL: `P0064_1002.sv.with_compound_het.tsv`, `P0064_1002.cnv.with_compound_het.tsv`. 