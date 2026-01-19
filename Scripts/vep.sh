#!/bin/bash

# Function: VEP annotation 
# Usage: sh vep.sh -i P001_110 -v P001_110.vcf.gz -o output/ -g GRCH38
# NewbornRisk using VEP v113, VarXOMics using VEP v114

# Default values
VEP_CACHE='/mnt/nas/Genomics/vep_cache'


# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory (required)
    -v, --vcf FILE                  Path for input VCF file (required)
    -g, --genome GENOME             Reference genome version, option: GRCH38, GRCH37 (default is GRCH38)
    --fork                          Threads for VEP
    -h, --help                      Display this help message

REQUIRED ARGUMENTS:
    - input-sample: Sample ID
    - output-directory: Output directory path
    - vcf: User uploaded vcf file
EOF
}

# Initialize variables
INPUT_SAMPLE=""
OUTPUT_DIR=""
VCF_FILE=""
GENOME_VERSION="GRCH38"
FORK="20"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-sample)
            INPUT_SAMPLE="$2"
            shift 2
            ;;
        -o|--output-directory)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -v|--vcf)
            VCF_FILE="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME_VERSION="$2"
            shift 2
            ;;
        --fork)
            FORK="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_SAMPLE" ]]; then
    echo "Error: --input-sample is required"
    usage
    exit 1
fi

if [[ -z "$VCF_FILE" ]]; then
    echo "Error: --input-sample is required"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-directory is required"
    usage
    exit 1
fi

if [[ -z "$GENOME_VERSION" ]]; then
    echo "Error: --genome is required"
    usage
    exit 1
fi

# Validate directories exist
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Error: Output directory does not exist: $OUTPUT_DIR"
    exit 1
fi


# Display configuration
echo "=== VEP Analysis Configuration ==="
echo "Input Sample: $INPUT_SAMPLE"
echo "Output Filename: $OUTPUT_DIR"
echo "VCF File: $VCF_FILE"
echo "Reference genome version is: $GENOME_VERSION"
echo "Fork for VEP $FORK"
echo "======================================"

printf "Start running vep annotation \ninput: $VCF_FILE  \n"
printf "vep cache: 113 \nvep version: 113 \n"


# Add 'chr' to chromosome if it is absence in the vcf file
BCFTOOLS="/mnt/nas/Genomics/miniconda3/envs/demultiplex/bin/bcftools"
CHR_FILE='/mnt/nas/Genomics/VarXOmics/tools/picard/chr_map.txt'
CHR_VCF="${VCF_FILE}_chr.vcf.gz"
$BCFTOOLS annotate --rename-chrs $CHR_FILE ${VCF_FILE} -Oz -o $CHR_VCF

# Define filenames
VEP_OUTPUT_VCF="$OUTPUT_DIR/${INPUT_SAMPLE}_vep_annotated.vcf.gz"
VEP_OUTPUT_BIALLELIC_VCF="$OUTPUT_DIR/${INPUT_SAMPLE}_vep_annotated_biallelic.vcf.gz"

#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate vep113, NewbornRisk using VEP v113
if [[ "$GENOME_VERSION" == "GRCH38" ]]; then
    echo  "Running VEP annotation on GRCH38. "
    vep --cache -dir_cache $VEP_CACHE \
        --offline \
        --cache_version 113 \
        --fork $FORK \
        --format vcf \
        --vcf \
        --dir_plugins $VEP_CACHE/VEP_plugins \
        -i $CHR_VCF \
        -o $VEP_OUTPUT_VCF \
        --force_overwrite \
        --compress_output bgzip \
        --assembly GRCh38 \
        --symbol --vcf --check_existing --variant_class \
        --sift b --polyphen b \
        --synonyms $VEP_CACHE/homo_sapiens_refseq/113_GRCh38/chr_synonyms.txt \
        --hgvs --refseq \
        --fasta $VEP_CACHE/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --canonical \
        --pick --pick_order mane_select,rank \
	    --af --af_gnomade --af_gnomadg --max_af \
        --custom $VEP_CACHE/vep_custom/clinvar_20250504.vcf.gz,ClinVar,vcf,exact,0,ID,CLNSIG,CLNDN,CLNHGVS,CLNSIGINCL,CLNVC,GENEINFO,CLNDISDB,CLNSIGCONF,CLNREVSTAT,CLNDNINCL \
        --plugin dbscSNV,$VEP_CACHE/dbscSNV1.1/dbscSNV1.1_GRCh38.txt.gz \
        --plugin REVEL,file=$VEP_CACHE/REVEL/new_tabbed_revel_grch38.tsv.gz \
        --plugin dbNSFP,$VEP_CACHE/dbNSFP4.7a/dbNSFP4.7a_grch38.gz,REVEL_score,REVEL_rankscore,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_rankscore,BayesDel_noAF_pred \
        --plugin SpliceAI,snv=$VEP_CACHE/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$VEP_CACHE/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin AlphaMissense,file=$VEP_CACHE/AlphaMissense/AlphaMissense_hg38.tsv.gz,cols=all \
        --verbose


elif [[ "$GENOME_VERSION" == "GRCH37" ]]; then
    echo "Reference genome is GRCh37. Liftover to GRCh38..."

    PICARD='/mnt/nas/Genomics/miniconda3/envs/vep114/share/picard-3.4.0-0/picard.jar'
    CHAIN_FILE="/mnt/nas/Genomics/VarXOmics/tools/picard/hg19ToHg38.over.chain"
    REF_FASTA="/mnt/nas/Genomics/VarXOmics/tools/picard/hg38.fa"

    # decompress VCF for liftover
    gunzip -c -f $CHR_VCF > "$OUTPUT_DIR/${INPUT_SAMPLE}_lifeover_tmp_chr.vcf" 

    # use GATK LiftoverVcf
    java -jar $PICARD LiftoverVcf \
        -I "$OUTPUT_DIR/${INPUT_SAMPLE}_lifeover_tmp_chr.vcf" \
        -O "$OUTPUT_DIR/${INPUT_SAMPLE}_lifted.vcf" \
        -CHAIN "$CHAIN_FILE" \
        -REJECT "$OUTPUT_DIR/${INPUT_SAMPLE}_rejected.vcf" \
        -R "$REF_FASTA"

    bgzip "$OUTPUT_DIR/${INPUT_SAMPLE}_lifted.vcf"
    mv "$OUTPUT_DIR/${INPUT_SAMPLE}_lifted.vcf.gz" $CHR_VCF
    #tabix -p vcf "$OUTPUT_DIR/${INPUT_SAMPLE}_PASS_lifted.vcf.gz"

    vep --cache -dir_cache $VEP_CACHE \
        --offline \
        --cache_version 113 \
        --fork $FORK \
        --format vcf \
        --vcf \
        --dir_plugins $VEP_CACHE/VEP_plugins \
        -i $CHR_VCF \
        -o $VEP_OUTPUT_VCF \
        --force_overwrite \
        --compress_output bgzip \
        --assembly GRCh38 \
        --symbol --vcf --check_existing --variant_class \
        --sift b --polyphen b \
        --synonyms $VEP_CACHE/homo_sapiens_refseq/113_GRCh38/chr_synonyms.txt \
        --hgvs --refseq \
        --fasta $VEP_CACHE/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --canonical \
        --pick --pick_order mane_select,rank \
	    --af --af_gnomade --af_gnomadg --max_af \
        --custom $VEP_CACHE/vep_custom/clinvar_20250504.vcf.gz,ClinVar,vcf,exact,0,ID,CLNSIG,CLNDN,CLNHGVS,CLNSIGINCL,CLNVC,GENEINFO,CLNDISDB,CLNSIGCONF,CLNREVSTAT,CLNDNINCL \
        --plugin dbscSNV,$VEP_CACHE/dbscSNV1.1/dbscSNV1.1_GRCh38.txt.gz \
        --plugin REVEL,file=$VEP_CACHE/REVEL/new_tabbed_revel_grch38.tsv.gz \
        --plugin dbNSFP,$VEP_CACHE/dbNSFP4.7a/dbNSFP4.7a_grch38.gz,REVEL_score,REVEL_rankscore,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_rankscore,BayesDel_noAF_pred \
        --plugin SpliceAI,snv=$VEP_CACHE/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$VEP_CACHE/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin AlphaMissense,file=$VEP_CACHE/AlphaMissense/AlphaMissense_hg38.tsv.gz,cols=all \
        --verbose

else
    echo "Error: Unsupported genome version '$GENOME_VERSION'. Only GRCH38 and GRCH37 are supported."
    exit 1
fi

# split into biallelic 
echo "Now spliting vep annotated vcf into biallelic file."
$BCFTOOLS norm -m - -Oz -o $VEP_OUTPUT_BIALLELIC_VCF $VEP_OUTPUT_VCF
rm $VEP_OUTPUT_VCF
mv $VEP_OUTPUT_BIALLELIC_VCF $VEP_OUTPUT_VCF 

echo -e "VEP finished at $(date).\n Python script will run now.  "
