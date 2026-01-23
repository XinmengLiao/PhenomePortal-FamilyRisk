### Commands 
```bash
# single - cnv/sv with snv
Scripts="/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts"
bash $Scripts/FamilyRisk_single_cnvsv.sh \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv/RapidWGS_001_C.cnv.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv \
	--genome GRCH38 \
	--snvindel-vcf /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_newborn/RapidWGS_001_C.hard-filtered.vcf.gz

# single - cnv/sv without snv
bash $Scripts/FamilyRisk_single_cnvsv.sh \
	-i rwgs1_kid \
	-v /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv/RapidWGS_001_C.cnv.vcf.gz \
	-o /mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv \
	--genome GRCH38
```

### UI design
Results are stored in `/mnt/nas/Genomics/Genome/FamilyRisk/examples20260119/single_cnvsv/Results`
It will be nice if the following contents could be included in the webpage:
1) Complete results of annotated CNV results: `RapidWGS_001_C.cnv.annotated.tsv`
2) Complete reullts of annotated CNV/SV with compound heterozygous variants identified with SNV/INDEL: `RapidWGS_001_C.cnv.annotated.CompoundHet.tsv`, `RapidWGS_001_C.sv.annotated.CompoundHet.tsv`