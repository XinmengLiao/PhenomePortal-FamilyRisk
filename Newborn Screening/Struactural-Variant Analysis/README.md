### Commands 
```bash
# single - cnv/sv with snv
bash Scripts/FamilyRisk_single_cnvsv.sh \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv/RapidWGS_001_C.cnv.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv \
	--genome GRCH38 \
	--snvindel-vcf /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_newborn/RapidWGS_001_C.hard-filtered.vcf.gz

# single - cnv/sv without snv
bash Scripts/FamilyRisk_single_cnvsv.sh \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv/RapidWGS_001_C.cnv.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv \
	--genome GRCH38
```
