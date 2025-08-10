## Input file 
- vcf file

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
- Detailed nodup table: `P0064_1203.txt`
- PGx `P0064_1203_PGx.txt`

## UI design
The UI could be similar to the original portal, selecting the variants, genes, and add comments and suggestions. 
