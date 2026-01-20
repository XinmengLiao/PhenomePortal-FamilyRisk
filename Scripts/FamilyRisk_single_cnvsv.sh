#!/bin/bash

# Script: CNV and SV annotation, including SNV/INDEL for determining the compound heterozygous state


# Default values
ANNOTSV='/mnt/nas/Genomics/VarXOmics/tools/AnnotSV/bin/AnnotSV'

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory path (required)
    -v, --vcf VCF_FILE              Input CNV/SV VCF file (required)
    --snvindel-vcf                  SNV/INDEL VCF file for compound heterozygous analysis (optional)
    -g, --genome GENOME             Reference Genome, default is GRCH38. 
                                    Options: GRCH37, GRCH38
    -h, --help                      Display this help message

REQUIRED ARGUMENTS:
    --input-sample: Sample ID
    --output-directory: Output directory path
    -v, --vcf: User uploaded vcf file
EOF
}

# Initialize variables
INPUT_SAMPLE=""
OUTPUT_DIR=""
VCF_FILE=""
GENOME='GRCH38'
SNV_INDEL_VCF=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-sample)
            INPUT_SAMPLE="$2"
            shift 2
            ;;
        -v|--vcf)
            VCF_FILE="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -o|--output-directory)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --snvindel-vcf)
            SNV_INDEL_VCF="$2"
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

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-directory is required"
    usage
    exit 1
fi

if [[ -z "$VCF_FILE" ]]; then
    echo "Error: -v is required"
    usage
    exit 1
fi


# Display configuration
echo "=== Exomiser Analysis Configuration ==="
echo "Input Sample: $INPUT_SAMPLE"
echo "Output Directory: $OUTPUT_DIR"
echo "Reference is $GENOME"
echo "Input VCF File: $VCF_FILE"
echo "======================================"


# Step 1: SV annotation
# Genome in AnnotSV is GRCh38 or GRCh37, not capital H

if [[ $GENOME == "GRCH38" ]]; then 
    AnnotSV_GENOME="GRCh38"
elif [[ $GENOME == "GRCH37" ]]; then
    AnnotSV_GENOME="GRCh37"
fi

mkdir -p "$OUTPUT_DIR"

if [[ -n "$SNV_INDEL_VCF" ]]; then
    echo "SNV/INDEL VCF provided: $SNV_INDEL_VCF"
    conda run -n vep114 $ANNOTSV \
        -SVinputFile "$VCF_FILE" \
        -outputDir "$OUTPUT_DIR" \
        -genomeBuild "$AnnotSV_GENOME" \
        -overwrite 1 -snvIndelPASS 1 \
        -candidateSnvIndelFiles "$SNV_INDEL_VCF"
    for file in "$OUTPUT_DIR"/*.annotated.tsv; do
        mv "$file" "${file/.annotated.tsv/.annotated.CompoundHet.tsv}"
    done
else
    echo "No SNV/INDEL VCF provided."
    conda run -n vep114 $ANNOTSV \
        -SVinputFile "$VCF_FILE" \
        -outputDir "$OUTPUT_DIR" \
        -genomeBuild "$AnnotSV_GENOME" \
        -overwrite 1 -snvIndelPASS 1
fi

rm -rf $OUTPUT_DIR/*.unannotated.tsv
rm -rf $OUTPUT_DIR/*AnnotSV_inputSV*.tsv
rm -rf $OUTPUT_DIR/*AnnotSV_inputSVfile.formatted.sorted*