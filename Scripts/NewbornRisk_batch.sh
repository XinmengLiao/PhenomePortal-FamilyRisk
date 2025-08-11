#!/bin/bash

# scripts for steamline the whole pipeline
SCRIPTS='/mnt/storage_pool/Genomics/Genome/NewbornRisk/Scripts'
#CONFIG_FILE="/mnt/storage_pool/Genomics/VarXOmics/examples/${1}/${1}_multivariant_config.yaml"


# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory path (required)
    -v, --vcf FILE                  Path for input VCF file (required)
    --sample-metadata               Metadata containing all Sample IDs (required)
                                    Format: [sampleID\tSex], Sex should be 'Male', 'Female', or blank
                                    i.e. P001_01    Female
    --genome						Reference genome, default is GRCH38
    								Options: GRCH37, GRCH38
    --only-pass                    	Only keep the PASS variant in vcf file. Default is yes. 
    								Options: yes, no. 
    --genedb                        The reference gene-disease list (required)
                                    Options: TR, babyseq, babydetect, babyscreen, guardian, earlycheck, ACMG, or user customized list
    -h, --help                      Display this help message
    --af-clinvar                    User defined allele frequency threshold for ClinVar variants, default is 1
    --af-precition                  User defined allele frequency threshold for predicted variants, default is 0.05
    --ada                           User defined ada threshold, default is 0.6
    --rf                            User defined rf threshold, default is 0.6
    --revel                         User defined revel threshold, default is 0.75
    --spliceai-al                   User defined SpliceAI acceptor loss threshold, default is 0.5
    --spliceai-ag                   User defined SpliceAI acceptor gain threshold, default is 0.5
    --spliceai-dl                   User defined SpliceAI donor loss threshold, default is 0.5
    --spliceai-dg                   User defined SpliceAI donor gain threshold, default is 0.5
    --bayesdel-addaf                User defined bayesdel (with AF) threshold, default is 0.0692655
    --bayesdel-noaf                 User defined bayesdel (no AF) threshold, default is 0.0570105
    --am-classification             User defined revel threshold, default is likely_pathogenic, ambigous
    --am-pathogenicity              User defined revel threshold, default is 0.564
    --clinvar                       User defined revel threshold, 
                                    default is Pathogenic, Likely_pathogenic, Uncertain_significance, Conflicting_classifications_of_pathogenicity
    --acmg-classification           User defined revel threshold, 
                                    default is Pathogenic, Likely_pathogenic, Uncertain_significance, Benign, Likely_benign

    --fork                          Threads for VEP, default is 20
    -t, --threads  THREADS          Threads for bcftools, default is 20 


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
GENOME="GRCH38"
ONLY_PASS="yes"
GENEDB=""
METADATA=""
FORK="20"
THREADS="20"

# Default thresholds
AF_CLINVAR="1"
AF_PRECITION="0.05"
ADA="0.6"
RF="0.6"
REVEL="0.75"
SPLICEAI_AL="0.5"
SPLICEAI_AG="0.5"
SPLICEAI_DL="0.5"
SPLICEAI_DG="0.5"
BAYESDEL_ADDAF="0.0692655"
BAYESDEL_NOAF="-0.0570105"
AM_CLASSIFICATION="likely_pathogenic,ambigous"
AM_PATHOGENICITY="0.564"
CLINVAR="Pathogenic,Likely_pathogenic,Uncertain_significance,Conflicting_classifications_of_pathogenicity"
ACMG_CLASSIFICATION="Pathogenic,Likely_pathogenic,Uncertain_significance,Benign,Likely_benign"


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
        --genome)
            GENOME="$2"
            shift 2
            ;;
        --genedb)
            GENEDB="$2"
            shift 2
            ;;
        --sample-metadata)
            METADATA="$2"
            shift 2
            ;;
        --only-pass)
            ONLY_PASS="$2"
            shift 2
            ;;            
        --af-clinvar)
            AF_CLINVAR="$2"
            shift 2
            ;;
        --af-precition)
            AF_PRECITION="$2"
            shift 2
            ;;
        --ada)
            ADA="$2"
            shift 2
            ;;
        --rf)
            RF="$2"
            shift 2
            ;;
        --revel)
            REVEL="$2"
            shift 2
            ;;
        --spliceai-al)
            SPLICEAI_AL="$2"
            shift 2
            ;;
        --spliceai-ag)
            SPLICEAI_AG="$2"
            shift 2
            ;;
        --spliceai-dl)
            SPLICEAI_DL="$2"
            shift 2
            ;;
        --spliceai-dg)
            SPLICEAI_DG="$2"
            shift 2
            ;;
        --bayesdel-addaf)
            BAYESDEL_ADDAF="$2"
            shift 2
            ;;
        --bayesdel-noaf)
            BAYESDEL_NOAF="$2"
            shift 2
            ;;
        --am-classification)
            AM_CLASSIFICATION="$2"
            shift 2
            ;;
        --am-pathogenicity)
            AM_PATHOGENICITY="$2"
            shift 2
            ;;
        --clinvar)
            CLINVAR="$2"
            shift 2
            ;;
        --acmg-classification)
            ACMG_CLASSIFICATION="$2"
            shift 2
            ;;
        --fork)
            FORK="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
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
    echo "Error: --vcf is required"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-directory is required"
    usage
    exit 1
fi

# Validate VCF file exists
if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: VCF file does not exist: $VCF_FILE"
    exit 1
fi

# Display configuration
echo "=== Sample Infomation ==="
echo "Input Sample: $INPUT_SAMPLE"
echo "VCF File: $VCF_FILE"
echo "Output Directory: $OUTPUT_DIR"
echo "Ped File: $PED"
echo "Reference Genome: $GENOME"
echo "Pass only: $ONLY_PASS"
echo "Gene-Disease Database: $GENEDB"
echo "File containing all sample IDs: $METADATA"
echo "VEP fork: $FORK"
echo "Bcftools threads: $THREADS"
echo "=== Threshold Parameters ==="
echo "AF ClinVar: $AF_CLINVAR"
echo "AF Precition: $AF_PRECITION"
echo "ADA: $ADA"
echo "RF: $RF"
echo "REVEL: $REVEL"
echo "SpliceAI AL: $SPLICEAI_AL"
echo "SpliceAI AG: $SPLICEAI_AG"
echo "SpliceAI DL: $SPLICEAI_DL"
echo "SpliceAI DG: $SPLICEAI_DG"
echo "BayesDel AddAF: $BAYESDEL_ADDAF"
echo "BayesDel NoAF: $BAYESDEL_NOAF"
echo "AM Classification: $AM_CLASSIFICATION"
echo "AM Pathogenicity: $AM_PATHOGENICITY"
echo "ClinVar: $CLINVAR"
echo "ACMG Classification: $ACMG_CLASSIFICATION"
echo "======================================"

mkdir -p $OUTPUT_DIR

### Step 1: Clean up the raw file 
echo '1. Remove missing ALT for input vcf. '
conda run -n vep bcftools view -e 'ALT = "."' $VCF_FILE -Oz -o $OUTPUT_DIR/${INPUT_SAMPLE}_rmmissingalt.vcf.gz --threads $THREADS
echo '2. Split into biallelic.'
conda run -n vep bcftools norm -m -both -Oz -o $OUTPUT_DIR/${INPUT_SAMPLE}_biallelic.vcf.gz $OUTPUT_DIR/${INPUT_SAMPLE}_rmmissingalt.vcf.gz --threads $THREADS
echo '3. Remove the duplicated variant'
conda run -n vep bcftools norm -d exact -Oz -o $OUTPUT_DIR/${INPUT_SAMPLE}_biallelic_nodup.vcf.gz $OUTPUT_DIR/${INPUT_SAMPLE}_biallelic.vcf.gz --threads $THREADS

# Step 2: VEP
echo "4. $(date) Run VEP for $INPUT_SAMPLE"
conda run -n vep bash $SCRIPTS/vep.sh \
    -v $OUTPUT_DIR/${INPUT_SAMPLE}_biallelic_nodup.vcf.gz \
    -o $OUTPUT_DIR \
    -i $INPUT_SAMPLE \
    -g $GENOME \
    --only_pass $ONLY_PASS \
    --fork $FORK

# Step 3: Python (Except GeneBe can not be run due to the too permerssive)
echo "5. $(date): Running Python for VEP result management..."

VEP_OUTPUT=$OUTPUT_DIR/${INPUT_SAMPLE}_vep_annotated.vcf.gz
mv ~/.netrc ~/.netrc.backup  # Remove /.netrc to remove the restriction from GeneBE
conda run -n vep python "$SCRIPTS/batch.py" \
    "$OUTPUT_DIR/${INPUT_SAMPLE}_vep_annotated.vcf.gz" \
    "$OUTPUT_DIR/${INPUT_SAMPLE}.txt" \
    "$METADATA" \
    "$AF_CLINVAR" \
    "$AF_PRECITION" \
    "$ADA" "$RF" \
    "$REVEL" \
    "$SPLICEAI_AL" \
    "$SPLICEAI_DG" \
    "$SPLICEAI_DL" \
    "$SPLICEAI_AG" \
    "$BAYESDEL_ADDAF" \
    "$BAYESDEL_NOAF" \
    "$AM_CLASSIFICATION" \
    "$AM_PATHOGENICITY" \
    "$CLINVAR" \
    "$ACMG_CLASSIFICATION" \
    "$GENEDB"
mv ~/.netrc.backup ~/.netrc

# Step 4: R for pedigree plot 
echo "6. Running R for creating pedigree plot"
mkdir -p $OUTPUT_DIR/Figures
conda run -n varxomics Rscript $SCRIPTS/Batch_Analysis.R $INPUT_SAMPLE $OUTPUT_DIR/${INPUT_SAMPLE}.txt $METADATA $OUTPUT_DIR $GENEDB 
# for biomRt package, the monogenic positive combined bar plot 
conda run -n varxomics2 Rscript $SCRIPTS/Batch_Analysis2.R $INPUT_SAMPLE $OUTPUT_DIR/${INPUT_SAMPLE}.txt $METADATA $OUTPUT_DIR $GENEDB 
