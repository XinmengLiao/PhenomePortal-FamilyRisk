#!/bin/bash

# scripts for steamline the whole pipeline
SCRIPTS='/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts'
DEFAULT_CARRIER_GENEDB='NBScreening'
DEFAULT_NEWBORN_GENEDB='ECS'

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory path (required)
    -v, --vcf FILE                  Path for input VCF file (required)
    --carrier                       Run Carrier Screening
    --newborn                       Run Newborn Risk Screening 
    --gender                        Gender ((required for PRS analysis)), 
                                    Options: Male, Female
    --genome						Reference genome, default is GRCH38
    								Options: GRCH37, GRCH38
    --only-pass                    	Only keep the PASS variant in vcf file. Default is yes. 
    								Options: yes, no. 
    --genedb                        The reference gene-disease list 
                                    Options: NBScreening, babyseq, babydetect, babyscreen, guardian, earlycheck, ACMG, ECS (Expanded Carrier Screening), or user customized list
                                    if --genedb is not provides, the gene-disese list will be set to NBScreening for newborn screening and expanded carrier screening list by default.
    --customized-genedb             User customized gene-disease list file path (required if genedb is customized) (optional)

For PRS analysis (optional):
    --run-prs                       Run PRS analysis or not (optional), default is no
                                    Options: yes, no
    --run-imputation                Run genotype imputation or not (optional). This argument is only for --run-prs.
                                    Options: yes, no, default is no    
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
    --input-sample: Sample ID
    --output-directory: Output directory path
    --vcf: User uploaded vcf file
    --carrier or --newborn: Choose either carrier screening or newborn risk screening
EOF
}

# Initialize variables
INPUT_SAMPLE=""
OUTPUT_DIR=""
VCF_FILE=""
GENOME="GRCH38"
ONLY_PASS="yes"
GENEDB=""
CUSTOMIZED_GENEDB=""
GENDER=""
RUNPRS="no"
RUNIMPUTATION="no"
FORK="20"
THREADS="20"
CARRIER=false
NEWBORN=false

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
        --customized-genedb)
            CUSTOMIZED_GENEDB="$2"
            shift 2
            ;;
        --gender)
            GENDER="$2"
            shift 2
            ;;
        --only-pass)
            ONLY_PASS="$2"
            shift 2
            ;;
        --carrier)
            CARRIER=true
            shift
            ;;
        --newborn)
            NEWBORN=true
            shift
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

# Check if either --carrier or --newborn is provided
if ! $CARRIER && ! $NEWBORN; then
    echo "Error: Either --carrier or --newborn must be provided"
    usage
    exit 1
fi

# Set default GENEDB based on flag if not provided by user
if [[ -z "$GENEDB" ]]; then
    if $CARRIER; then
        GENEDB="$DEFAULT_CARRIER_GENEDB"
    elif $NEWBORN; then
        GENEDB="$DEFAULT_NEWBORN_GENEDB"
    fi
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
echo "Customized Gene-Disease Database: $CUSTOMIZED_GENEDB"
echo "Sample Gender: $GENDER"
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

# File names definiation
INPUT_VCF_RMMISSINGALT="${INPUT_SAMPLE}_rmmissingalt.vcf.gz"
INPUT_VCF_BIALLELIC="${INPUT_SAMPLE}_biallelic.vcf.gz"
INPUT_VCF_BIALLELIC_NODUP="${INPUT_SAMPLE}_biallelic_nodup.vcf.gz"
INPUT_VCF_BIALLELIC_NODUP_PASS="${INPUT_SAMPLE}_biallelic_nodup_pass.vcf.gz"
INPUT_VCF_ANNOTATED="${INPUT_SAMPLE}_vep_annotated.vcf.gz"

# ### Step 1: Clean up the raw file 
# echo '1. Remove missing ALT for input vcf. '
# conda run -n vep bcftools view -e 'ALT = "."' $VCF_FILE -Oz -o $OUTPUT_DIR/${INPUT_VCF_RMMISSINGALT} --threads $THREADS
# echo '2. Split into biallelic.'
# conda run -n vep bcftools norm -m -both -Oz -o $OUTPUT_DIR/${INPUT_VCF_BIALLELIC} $OUTPUT_DIR/${INPUT_VCF_RMMISSINGALT} --threads $THREADS
# echo '3. Remove the duplicated variant'
# conda run -n vep bcftools norm -d exact -Oz -o $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP} $OUTPUT_DIR/${INPUT_VCF_BIALLELIC} --threads $THREADS

# ### Step 2: Filtered sequencing quality PASS
# echo "3. $(date) Filtered sequencing quality PASS"
# if [[ "$ONLY_PASS" == "yes" ]]; then
#     zgrep -E "^#|PASS" $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP} | bgzip > $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS}
# else
#     mv $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP} $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS}
# fi

# ### Step 2.5 Remove intermediate files
# rm $OUTPUT_DIR/${INPUT_VCF_RMMISSINGALT} $OUTPUT_DIR/${INPUT_VCF_BIALLELIC} $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP}
# echo 'Finished cleaning up the input vcf file.'

# Step 3: VEP
echo "4. $(date) Run VEP for $INPUT_SAMPLE"
conda run -n vep bash $SCRIPTS/vep.sh \
    -v $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS} \
    -o $OUTPUT_DIR \
    -i $INPUT_SAMPLE \
    -g $GENOME \
    --fork $FORK

# # Step 4: Python (Except GeneBe can not be run due to the too permerssive)
# echo "5. $(date): Running Python for VEP result management..."

# conda run -n vep python "$SCRIPTS/single.py" \
# 	"$OUTPUT_DIR/${INPUT_VCF_ANNOTATED}" \
#     "$OUTPUT_DIR/${INPUT_SAMPLE}.txt" \
#     "$GENDER" \
#     "$AF_CLINVAR" \
#     "$AF_PRECITION" \
#     "$ADA" "$RF" \
#     "$REVEL" \
#     "$SPLICEAI_AL" \
#     "$SPLICEAI_DG" \
#     "$SPLICEAI_DL" \
#     "$SPLICEAI_AG" \
#     "$BAYESDEL_ADDAF" \
#     "$BAYESDEL_NOAF" \
#     "$AM_CLASSIFICATION" \
#     "$AM_PATHOGENICITY" \
#     "$CLINVAR" \
#     "$ACMG_CLASSIFICATION" \
#     "$GENEDB" \
#     "$CUSTOMIZED_GENEDB"

# # Step 5: PGx by PharmCat
# pharmcat="/mnt/nas/Genomics/Genome/NewbornRisk/tools/pharmcat/pharmcat-3.1.1-all.jar"
# pharmcat_preprocessor="/mnt/nas/Genomics/Genome/NewbornRisk/tools/pharmcat/preprocessor/pharmcat_vcf_preprocessor"
# preprocessor_ref="/mnt/nas/Genomics/Genome/NewbornRisk/tools/pharmcat/reference.fna.bgz"
# preprocessor_position="/mnt/nas/Genomics/Genome/NewbornRisk/tools/pharmcat/pharmcat_positions_3.1.1.vcf.bgz"

# PHARMCAT_PREPROCESSED_VCF="${INPUT_SAMPLE}_biallelic_nodup_pass.preprocessed.vcf.bgz"

# mkdir -p ${OUTPUT_DIR}/PGx

# # normalized by pharmcat preprocessor
# conda run -n vep114 $pharmcat_preprocessor \
#     -vcf $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS} \
#     -refFna $preprocessor_ref \
#     -refVcf $preprocessor_position

# # pharmcat step A  
# java -jar $pharmcat \
#     -matcher -vcf $PHARMCAT_PREPROCESSED_VCF \
#     -phenotyper -o ${OUTPUT_DIR}/PGx/ \
#     -research cyp2d6 -v

# # pharmcat step B
# java -jar $pharmcat \
#     -reporter -ri ${OUTPUT_DIR}/PGx/*.phenotype.json \
#     -o ${OUTPUT_DIR}/PGx/ -reporterJson -reporterHtml

# # Step 6: PRS analysis (optional)
# if [[ "$RUNPRS" == "yes" ]]; then
#     echo "$(date): Running PRS analysis for $INPUT_SAMPLE ..."
#     bash $SCRIPTS/NewbornRisk_PRS_Single.sh \
#         -i $INPUT_SAMPLE \
#         -v $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS} \    
#         -o $OUTPUT_DIR \                
#         --genome $GENOME \
#         --gender $GENDER \
#         --run-imputation $RUNIMPUTATION \
#         $ID_ARG
#         --only-pass $ONLY_PASS \
#         -t $THREADS
# else
#     echo "$(date): Skipping PRS analysis for $INPUT_SAMPLE ..."
# fi  