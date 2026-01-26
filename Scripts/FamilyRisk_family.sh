#!/bin/bash

# scripts for steamline the whole pipeline
SCRIPTS='/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts'

# Exit immediately if any command fails
set -e


# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory path (required)
    -v, --vcf FILE                  Path for input VCF file for family-merged vcf file
    --sample-list LIST              Line separated list (tsv/txt) of multiple family vcf files, which would be merged by FamilyRisk
    --carrier                       Run Carrier Screening
    --newborn                       Run Newborn Risk Screening 
    --ped                           Ped file containing family relationships (required)
    --genome						Reference genome, default is GRCH38
    								Options: GRCH37, GRCH38
    --only-pass                    	Only keep the PASS variant in vcf file. Default is yes. 
    								Options: yes, no. 
    --genedb                        The reference gene-disease list (required)
                                    Options: TR, babyseq, babydetect, babyscreen, guardian, earlycheck, ACMG, or user customized list
    --customized-genedb             The user customized gene-disease list file path

    For PGx and PRS analysis (optional):
    --run-pgx                       Run PGx analysis or not (optional), default is no
                                    Options: yes, no
    --pgsid						    A comma separated list of PGS score IDs, e.g. PGS000802
    --pgpid						    A comma separated list of PGS Catalog publications, e.g. PGP000001
    --efoid						    A comma separated list of PGS Catalog EFO traits, e.g. EFO_0004214
    --run-prs                       Whether to run PRS analysis, default is no. 
                                    Options: yes, no.
    --run-imputation                Whether to run imputation step, default is no. Yes is recommended for vcf files without 0/0 genotypes for PRS analysis.
                                    Options: yes, no.
    
    Variants Filteration Parameters
    --only-clinvar                  Only keep ClinVar reported variants for analysis, default is no.
                                    Options: yes, no.   
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
    
    Running Parameters
    --fork                          Threads for VEP, default is 20
    -t, --threads  THREADS          Threads for bcftools, default is 20 
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
PED=""
GENOME="GRCH38"
ONLY_PASS="yes"
GENEDB=""
CUSTOMIZED_GENEDB=""
FORK=""
THREADS=""
CARRIER=false
NEWBORN=false

# PRS parameters
RUNPRS="no"
RUNIMPUTATION="no"
PGSID=""
PGPID=""
EFOID=""

# PGx parameters
RUNPGX="no"

# Default thresholds
ONLY_CLINVAR="no"
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
        --sample-list)
            SAMPLE_LIST="$2"
            shift 2
            ;;
        --ped)
            PED="$2"
            shift 2
            ;;
        --genome)
            GENOME="$2"
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
        --only-pass)
            ONLY_PASS="$2"
            shift 2
            ;; 
        --only-clinvar)
            ONLY_CLINVAR="$2"
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
        --run-pgx)
            RUNPGX="$2"
            shift 2
            ;;
        --run-prs)
            RUNPRS="$2"
            shift 2
            ;;
        --run-imputation)
            RUNIMPUTATION="$2"
            shift 2
            ;;
        --pgsid)
            PGSID="$2"
            shift 2
            ;;
        --pgpid)
            PGPID="$2"
            shift 2
            ;;
        --efoid)
            EFOID="$2"
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

if [[ -z "$VCF_FILE" && -z "$SAMPLE_LIST" ]]; then
    echo "Error: either --vcf or --sample-list is required"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-directory is required"
    usage
    exit 1
fi

if [[ -z "$PED" ]]; then
    echo "Error: --ped is required"
    usage
    exit 1
fi

# Check if either --carrier or --newborn is provided
if ! $CARRIER && ! $NEWBORN; then
    echo "Error: Either --carrier or --newborn must be provided"
    usage
    exit 1
fi

# Set FUNCTION_TYPE based on flag
if $CARRIER; then
    FUNCTION_TYPE="carrier"
elif $NEWBORN; then
    FUNCTION_TYPE="newborn"
fi

# Set default GENEDB based on flag if not provided by user
if [[ -z "$GENEDB" ]]; then
    if $CARRIER; then
        GENEDB="$DEFAULT_CARRIER_GENEDB"
    elif $NEWBORN; then
        GENEDB="$DEFAULT_NEWBORN_GENEDB"
    fi
fi

# Display configuration
echo "=== Sample Infomation ==="
echo "Input Sample: $INPUT_SAMPLE"
echo "VCF File: $VCF_FILE"
echo "Sample List: $SAMPLE_LIST"
echo "Output Directory: $OUTPUT_DIR"
echo "Ped File: $PED"
echo "Reference Genome: $GENOME"
echo "Pass only: $ONLY_PASS"
echo "Gene-Disease Database: $GENEDB"
echo "Customized Gene-Disease Database: $CUSTOMIZED_GENEDB"
echo "VEP fork: $FORK"
echo "Bcftools threads: $THREADS"

echo "=== PRS and PGx Parameters ==="
echo "Run PRS analysis: $RUNPRS"
echo "Run Imputation: $RUNIMPUTATION"
echo "PGS IDs: $PGSID"
echo "PGS Publications: $PGPID"
echo "PGS EFO Traits: $EFOID"
echo "Run PGx analysis: $RUNPGX"

echo "=== Threshold Parameters ==="
echo "Only ClinVar: $ONLY_CLINVAR"
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

mkdir -p $OUTPUT_DIR/Results

## File names definiation
INPUT_VCF_RMMISSINGALT="${INPUT_SAMPLE}_rmmissingalt.vcf.gz"
INPUT_VCF_BIALLELIC="${INPUT_SAMPLE}_biallelic.vcf.gz"
INPUT_VCF_BIALLELIC_NODUP="${INPUT_SAMPLE}_biallelic_nodup.vcf.gz"
INPUT_VCF_BIALLELIC_NODUP_PASS="${INPUT_SAMPLE}_biallelic_nodup_pass.vcf.gz"
INPUT_VCF_ANNOTATED="${INPUT_SAMPLE}_vep_annotated.vcf.gz"
VCF_LIST="${OUTPUT_DIR}/${INPUT_SAMPLE}_vcf_list.txt"

### Step 0: If multiple sample vcf files are provided, merge them into one family vcf file
if [[ -f "$VCF_FILE" && -s "$VCF_FILE" ]]; then
    echo "User uploaded merged family VCF file: $VCF_FILE"
elif [[ -f "$SAMPLE_LIST" && -s "$SAMPLE_LIST" ]]; then
    echo "Merging VCF files from $SAMPLE_LIST..."
    
    while IFS= read -r vcf_path; do
        # skip blank and commented lines
        [[ -z "$vcf_path" || "$vcf_path" =~ ^# ]] && continue
        
        # Check if VCF file exists
        if [[ -f "$vcf_path" ]]; then
            echo "Found: $vcf_path"
            echo "$vcf_path" >> "$VCF_LIST"
            tabix -f -p vcf "$vcf_path"
        else
            echo "Warning: VCF file not found: $vcf_path, skip in merging step. "
        fi
    done < "$SAMPLE_LIST"
    
    # Merge all VCF files
    MERGED_VCF="${OUTPUT_DIR}/${INPUT_SAMPLE}_merged.vcf.gz"
    conda run -n vep bcftools merge -l "$VCF_LIST" -Oz -o "$MERGED_VCF" --threads "$THREADS"
    VCF_FILE="$MERGED_VCF"

    rm "$VCF_LIST"

else
    echo "Error: neither VCF_FILE nor SAMPLE_LIST is valid"
    exit 1
fi

### Step 1: Clean up the raw file 
echo '1. Remove missing ALT for input vcf. '
conda run -n vep bcftools view -e 'ALT = "."' $VCF_FILE -Oz -o $OUTPUT_DIR/${INPUT_VCF_RMMISSINGALT} --threads $THREADS
echo '2. Split into biallelic.'
conda run -n vep bcftools norm -m -both -Oz -o $OUTPUT_DIR/${INPUT_VCF_BIALLELIC} $OUTPUT_DIR/${INPUT_VCF_RMMISSINGALT} --threads $THREADS
echo '3. Remove the duplicated variant'
conda run -n vep bcftools norm -d exact -Oz -o $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP} $OUTPUT_DIR/${INPUT_VCF_BIALLELIC} --threads $THREADS

### Step 2: Filtered sequencing quality PASS
echo "3. $(date) Filtered sequencing quality PASS"
if [[ "$ONLY_PASS" == "yes" ]]; then
    zgrep -E "^#|PASS" $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP} | bgzip > $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS}
else
    mv $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP} $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS}
fi

### Step 2.5 Remove intermediate files
rm $OUTPUT_DIR/${INPUT_VCF_RMMISSINGALT} $OUTPUT_DIR/${INPUT_VCF_BIALLELIC} $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP}
echo 'Finished cleaning up the input vcf file.'

### Step 3: VEP
echo "4. $(date) Run VEP for $INPUT_SAMPLE"
conda run -n vep bash $SCRIPTS/vep.sh \
    -v $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS} \
    -o $OUTPUT_DIR \
    -i $INPUT_SAMPLE \
    -g $GENOME \
    --fork $FORK

### Step 4: Python 
# The log file is in $OUTPUT_DIR/family_debug.log
echo "5. $(date): Running Python for VEP result management..."

conda run -n vep python "$SCRIPTS/family.py" \
    "$OUTPUT_DIR/${INPUT_VCF_ANNOTATED}" \
    "$OUTPUT_DIR/${INPUT_SAMPLE}.txt" \
    "$PED" \
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
    "$GENEDB" \
    "$CUSTOMIZED_GENEDB" \
    "$FUNCTION_TYPE" \
    "$ONLY_CLINVAR"


### Step 5: R for pedigree plot 
echo "6. Running R for managing data and creating pedigree plot. "
if [[ $FUNCTION_TYPE == "carrier" ]]; then
    echo "Running Carrier Screening R script. "
    conda run -n varxomics Rscript $SCRIPTS/Carrier_Family20260121.R $INPUT_SAMPLE $OUTPUT_DIR/${INPUT_SAMPLE}.txt $PED $OUTPUT_DIR $GENEDB
elif [[ $FUNCTION_TYPE == "newborn" ]]; then
    echo "Running Newborn Risk Screening R script. "
    conda run -n varxomics Rscript $SCRIPTS/Newborn_Family20260121.R $INPUT_SAMPLE $OUTPUT_DIR/${INPUT_SAMPLE}.txt $PED $OUTPUT_DIR $GENEDB
fi
mv $OUTPUT_DIR/${INPUT_SAMPLE}.txt $OUTPUT_DIR/Results/

### Step 6: PGx by PharmCat
if [[ "$RUNPGX" == "no" ]]; then
    echo "Skipping PGx analysis as --run-pgx is not set to 'yes'."
    exit 0
else
    pharmcat="/mnt/nas/Genomics/Genome/FamilyRisk/tools/pharmcat/pharmcat-3.1.1-all.jar"
    pharmcat_preprocessor="/mnt/nas/Genomics/Genome/FamilyRisk/tools/pharmcat/preprocessor/pharmcat_vcf_preprocessor"
    preprocessor_ref="/mnt/nas/Genomics/Genome/FamilyRisk/tools/pharmcat/reference.fna.bgz"
    preprocessor_position="/mnt/nas/Genomics/Genome/FamilyRisk/tools/pharmcat/pharmcat_positions_3.1.1.vcf.bgz"

    PHARMCAT_PREPROCESSED_VCF="${OUTPUT_DIR}/${INPUT_SAMPLE}_biallelic_nodup_pass.preprocessed.vcf.bgz"
    echo $PHARMCAT_PREPROCESSED_VCF

    mkdir -p ${OUTPUT_DIR}/PGx
    mkdir -p ${OUTPUT_DIR}/Results/PGx_Reports

    # normalized by pharmcat preprocessor
    conda run -n vep114 $pharmcat_preprocessor \
        -vcf $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS} \
        -refFna $preprocessor_ref \
        -refVcf $preprocessor_position

    # pharmcat step A  
    java -jar $pharmcat \
        -matcher -vcf $PHARMCAT_PREPROCESSED_VCF \
        -phenotyper -o ${OUTPUT_DIR}/PGx/ \
        -research cyp2d6 -v

    # pharmcat step B
    for file in ${OUTPUT_DIR}/PGx/*.phenotype.json; do
        echo $file
        java -jar $pharmcat \
            -reporter -ri $file \
            -o "${OUTPUT_DIR}/Results/PGx_Reports" -reporterJson -reporterHtml -v
    done

    # Remove intermediate files
    rm -f "${PHARMCAT_PREPROCESSED_VCF}"*
    rm -f "${OUTPUT_DIR}/PGx"/*missing_pgx_var.vcf
    rm -f "${OUTPUT_DIR}"/*missing_pgx_var.vcf

fi

### Step 7: Run PRS analysis if required
if [[ "$RUNPRS" == "yes" ]]; then
    
    if [[ -z "$RUNIMPUTATION" ]]; then
        echo "Error: --run-imputation must be set when --run-prs is 'yes'"
        exit 1
    fi

    if [[ -z "$GENDER" ]]; then
        echo "Error: --gender is required for PRS analysis"
        exit 1
    fi

    echo "7. Running PRS analysis for $INPUT_SAMPLE"

    # make PRS psam file 
    PSAM="$OUTPUT_DIR/${INPUT_SAMPLE}.psam"
    rm -f "$PSAM"
    printf "#IID\tSEX\n" >> "$PSAM"
    # extract individual ID (column 2) and sex (column 5) for the given sample
    awk '{print $1"\t"$2}' $METADATA >> $PSAM

    mkdir -p $OUTPUT_DIR/PRS

    echo "Now running PGS score for: $INPUT_VCF_BIALLELIC_NODUP_PASS"

    if [[ -n "$PGSID" ]]; then
    echo "Use PGSID: $PGSID. "
    conda run -n pgsc bash $SCRIPTS/FamilyRisk_PRS_Batch.sh \
        -i $INPUT_SAMPLE \
        -o $OUTPUT_DIR/PRS \
        -v $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP_PASS} \
        --metadata $PSAM \
        --genome $GENOME \
        --only-pass no \
        --run-imputation $RUNIMPUTATION \
        --pgsid $PGSID \
        -t $THREADS || { echo "PRS analysis failed"; exit 1; }

    elif [[ -n "$PGPID" ]]; then
        echo "Use PGPID: $PGPID. "
        conda run -n pgsc bash $SCRIPTS/FamilyRisk_PRS_Batch.sh \
            -i $INPUT_SAMPLE \
            -o $OUTPUT_DIR/PRS \
            -v $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP} \
            --metadata $PSAM \
            --genome $GENOME \
            --only-pass no \
            --run-imputation $RUNIMPUTATION \
            --pgpid $PGPID \
            -t $THREADS || { echo "PRS analysis failed"; exit 1; }

    elif [[ -n "$EFOID" ]]; then
        echo "Use EFOID: $EFOID. "
        conda run -n pgsc bash $SCRIPTS/FamilyRisk_PRS_Batch.sh \
            -i $INPUT_SAMPLE \
            -o $OUTPUT_DIR/PRS \
            -v $OUTPUT_DIR/${INPUT_VCF_BIALLELIC_NODUP} \
            --metadata $PSAM \
            --genome $GENOME \
            --only-pass no \
            --run-imputation $RUNIMPUTATION \
            --efoid $EFOID \
            -t $THREADS || { echo "PRS analysis failed"; exit 1; }

    fi

    # move Reports to the Results folder
    mv $OUTPUT_DIR/PRS/results/$INPUT_SAMPLE/score/ $OUTPUT_DIR/Results/PGS_Scores

fi

## Step 8: R for PGS density plot 
if [[ "$RUNPRS" == "yes" ]]; then
    echo "Running R for PGS density plots."
    gunzip -d -k $OUTPUT_DIR/Results/PGS_Scores/*popsimilarity.txt.gz
    gunzip -d -k $OUTPUT_DIR/Results/PGS_Scores/*pgs.txt.gz
    conda run -n varxomics2 Rscript $SCRIPTS/PGS_DensityPlot20260121.R \
        $INPUT_SAMPLE \
        $OUTPUT_DIR/Results/PGS_Scores/*pgs.txt.gz \
        $OUTPUT_DIR/Results/PGS_Scores/*popsimilarity.txt \
        $OUTPUT_DIR family
fi