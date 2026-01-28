#!/bin/bash

# This script is for running PRS for merged VCF files, including cohort and family-based analysis.
#set -e
#set -o pipefail  # Exit on pipe failure too

# scripts for steamline the whole pipeline
beagle="/mnt/nas/Genomics/Genome/FamilyRisk/tools/beagle5.jar"
nextflow="/mnt/nas/Genomics/Genome/FamilyRisk/tools/nextflow"
assembly38='/mnt/nas/Genomics/vep_cache/vep114_assembly/Homo_sapiens.GRCh38.dna.toplevel.fa.gz'
assembly37='/mnt/nas/Genomics/vep_cache/vep114_assembly/Homo_sapiens.GRCh37.dna.toplevel.fa.gz'
map_dir38='/mnt/nas/Genomics/Genome/FamilyRisk/Datasets/Beagle/plink.GRCh38/no_chr_in_chrom_field'
map_dir37='/mnt/nas/Genomics/Genome/FamilyRisk/Datasets/Beagle/plink.GRCh37/no_chr_in_chrom_field'
ref_dir38='/mnt/nas/Genomics/1KGenome/GRCh38/no_chr_in_chrom_field'
ref_dir37='/mnt/nas/Genomics/1KGenome/GRCh37-Phase3/no_chr_in_chrom_field'
remove_chr_file='/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts/remove_chrs.txt'
contig_header='/mnt/nas/Genomics/Genome/FamilyRisk/PhenomePortal-FamilyRisk/Scripts/contigs.txt'
# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    -i, --input-sample SAMPLE       Input sample name (required)
    -o, --output-directory DIR      Output directory path (required)
    -v, --vcf FILE                  Path for input VCF file (required)
    --pgsid						    A comma separated list of PGS score IDs, e.g. PGS000802
    --pgpid						    A comma separated list of PGS Catalog publications, e.g. PGP000001
    --efoid						    A comma separated list of PGS Catalog EFO traits, e.g. EFO_0004214
    --metadata                      Sample metadata in .txt or .tsv format (required), formatting with #IID,SEX for all newborns in the vcf files. 
                                    #IID should be the same as sample ID in vcf file.
                                    SEX should be 1 for male and 2 for female.    
    --genome						Reference genome, default is GRCH38
    								Options: GRCH37, GRCH38
    --run-imputation                Whether to run imputation step, default is yes. 
                                    Options: yes, no. 
                                    It is recommended to run imputation for WGS vcf files with 0/0 genotypes. 
    -h, --help                      Display this help message
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
SAMPLE_META=""
GENOME="GRCH38"
THREADS="20"  # Default value
PGSID=""
PGPID=""
EFOID=""
RUNIMPUTATION="yes"

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
        --metadata)
            SAMPLE_META="$2"
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
        --run-imputation)
            RUNIMPUTATION="$2"  
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


count=0
[[ -n "$PGSID" ]] && ((count++))
[[ -n "$EFOID" ]] && ((count++))
[[ -n "$PGPID" ]] && ((count++))

if [[ $count -eq 0 ]]; then
    echo "Error: one of --pgsid, --pgpid, or --efoid is required"
    usage
    exit 1
elif [[ $count -gt 1 ]]; then
    echo "Error: only ONE of --pgsid, --pgpid, or --efoid can be specified"
    usage
    exit 1
fi


# Display configuration
echo "=== Sample Infomation ==="
echo "Input Sample ID: $INPUT_SAMPLE"
echo "VCF File: $VCF_FILE"
echo "Output Directory: $OUTPUT_DIR"
echo "Reference Genome: $GENOME"
echo "Sample Metadata: $SAMPLE_META"
echo "Run Imputation: $RUNIMPUTATION"
echo "PGS ID: $PGSID"
echo "PGP ID: $PGPID"
echo "EFO ID: $EFOID"
echo "========================="

mkdir -p $OUTPUT_DIR

#remove the underscore in INPUT_SAMPLE
newname=$(echo "$INPUT_SAMPLE" | sed 's/_/-/g')
INPUT_SAMPLE=$newname

# set reference genome
if [[ "$GENOME" == "GRCH38" ]]; then
    assembly=$assembly38
    map_dir=$map_dir38
    ref_dir=$ref_dir38
    target_build="GRCh38"
elif [[ "$GENOME" == "GRCH37" ]]; then
    assembly=$assembly37
    map_dir=$map_dir37
    ref_dir=$ref_dir37
    target_build="GRCh37"
else
    echo "Error: Unsupported genome build $GENOME. Only GRCH38 is supported."
    exit 1
fi

# # Step 1: remove the 'chr' in users VCF file
echo "Removing 'chr' in $VCF_FILE at $(date)." 
bcftools annotate --rename-chrs $remove_chr_file $VCF_FILE -Oz -o $OUTPUT_DIR/${INPUT_SAMPLE}-renamed.vcf.gz --threads $THREADS || { echo "ERROR: bcftools annotate failed"; exit 1; }
tabix -p vcf $OUTPUT_DIR/${INPUT_SAMPLE}-renamed.vcf.gz -f

# Step 2: separate into 1-22 chromosomes for imputation (If imputation is required)
if [[ "$RUNIMPUTATION" != "yes" ]]; then
    echo "Skipping imputation step as per user request."
    mv $OUTPUT_DIR/${INPUT_SAMPLE}-renamed.vcf.gz $OUTPUT_DIR/${INPUT_SAMPLE}-imputed.allchr.vcf.gz

    printf "sampleset,path_prefix,chrom,format\n" > $OUTPUT_DIR/chr${i}/pgsc-input.csv
    line=$(echo "${INPUT_SAMPLE},$OUTPUT_DIR/${INPUT_SAMPLE}-imputed.allchr,,pfile" | sed 's/_/-/g')
    printf "${line}\n" >> $OUTPUT_DIR/pgsc-input.csv
else
    echo "Running imputation step from $(date)."

    # Read sample metadata and convert into psam file
    echo "Writing the psam file at $(date)."
    mv $SAMPLE_META $OUTPUT_DIR/${INPUT_SAMPLE}-imputed.allchr.psam
    printf "sampleset,path_prefix,chrom,format\n" > $OUTPUT_DIR/${INPUT_SAMPLE}-pgsc-input.csv

    for i in $(seq 1 22); do 
    #for i in 1; do
        SUB_OUTPUT_DIR=$OUTPUT_DIR/chromosome-file/chr${i}
        mkdir -p $SUB_OUTPUT_DIR
        
        # define filenames
        SPLIT_VCF=$SUB_OUTPUT_DIR/${INPUT_SAMPLE}_chr${i}.vcf.gz
        NORM_VCF=$SUB_OUTPUT_DIR/${INPUT_SAMPLE}_chr${i}_norm.vcf.gz
        IMPUTE_VCF=$SUB_OUTPUT_DIR/${INPUT_SAMPLE}_chr${i}_norm_impute.vcf.gz
        IMPUTE_FILTERED_VCF=$SUB_OUTPUT_DIR/${INPUT_SAMPLE}_chr${i}_norm_impute_filtered.vcf.gz
        REHEADER_VCF=$SUB_OUTPUT_DIR/${INPUT_SAMPLE}_chr${i}_norm_impute_filtered_reheader.vcf.gz
        RMDUP_VCF=$SUB_OUTPUT_DIR/${INPUT_SAMPLE}_chr${i}_rmdup.vcf.gz

        # split into 1-22
        echo "Split into chr${i} at $(date)."
        bcftools view -r ${i} $OUTPUT_DIR/${INPUT_SAMPLE}-renamed.vcf.gz \
            -Oz -o $SPLIT_VCF --threads $THREADS
        tabix -p vcf $SPLIT_VCF -f
    
        # normalize with assembly
        echo "Normalized $SPLIT_VCF to $GENOME at $(date)."
        bcftools norm -f $assembly -m -both \
            $SPLIT_VCF \
            -Oz -o $NORM_VCF --threads $THREADS
    
        # impute by beagle
        echo "$(date): Imputing $NORM_VCF."
        java -jar $beagle \
            gt=$NORM_VCF \
            ref=$ref_dir/ALL.chr${i}.vcf.gz \
            out=$SUB_OUTPUT_DIR/${INPUT_SAMPLE}_chr${i}_norm_impute \
            map=$map_dir/plink.chr${i}.map \
            nthreads=$THREADS
        
        # filtered out low-quality imputed variants
        echo "$(date): Filtering out low-quality variants for $IMPUTE_VCF"
        bcftools view -i '(GT="0|0") || (INFO/DR2>0.3)' \
            $IMPUTE_VCF \
            -Oz -o $IMPUTE_FILTERED_VCF \
            --threads $THREADS

        rm -rf $SUB_OUTPUT_DIR/header.txt
        bcftools view -h $IMPUTE_FILTERED_VCF | head -n -1 >> $SUB_OUTPUT_DIR/header.txt
        cat $contig_header >> $SUB_OUTPUT_DIR/header.txt
        bcftools view -h $IMPUTE_FILTERED_VCF | tail -n 1 >> $SUB_OUTPUT_DIR/header.txt
        
        bcftools reheader -h $SUB_OUTPUT_DIR/header.txt $IMPUTE_FILTERED_VCF > $REHEADER_VCF
        tabix -p vcf $REHEADER_VCF -f
        
        bcftools annotate $REHEADER_VCF \
            --set-id='%CHROM:%POS:%REF:%ALT' --threads $THREADS \
            -Oz -o $RMDUP_VCF

        # Convert each chromosome's vcf into the pfile for pgsc_calc
        echo "Converting chr${i} into pfile at $(date)."
        plink2 --vcf $RMDUP_VCF \
            --allow-extra-chr --threads $THREADS \
            --chr 1-22 \
            --make-pgen \
            --out $SUB_OUTPUT_DIR/${INPUT_SAMPLE}-chr${i}-imputed.allchr \
            --psam $OUTPUT_DIR/${INPUT_SAMPLE}-imputed.allchr.psam \
            --new-id-max-allele-len 10000 # unlimit the maximum ID length

        #Write input.csv
        echo "Writing pgsc-input.csv at $(date)."
        printf "sampleset,path_prefix,chrom,format\n" > $SUB_OUTPUT_DIR/pgsc-input.csv
        line=$(echo "${INPUT_SAMPLE},$SUB_OUTPUT_DIR/${INPUT_SAMPLE}-chr${i}-imputed.allchr,${i},pfile" | sed 's/_/-/g')
        printf "${line}\n" >> $SUB_OUTPUT_DIR/pgsc-input.csv

        echo "Combining all chromosomes' pgsc-input.csv."
        tail -n 1 $SUB_OUTPUT_DIR/pgsc-input.csv >> $OUTPUT_DIR/${INPUT_SAMPLE}-pgsc-input.csv

    done

fi

# Step 3: Run PRS with pgsc_calc
echo "Running pgsc_calc at $(date)."

cd $OUTPUT_DIR

if [[ -n "$PGSID" ]]; then
    echo "Use PGSID: $PGSID. "
    $nextflow run pgscatalog/pgsc_calc \
        --input $OUTPUT_DIR/${INPUT_SAMPLE}-pgsc-input.csv \
        --pgs_id $PGSID \
        --target_build $target_build \
        --min_overlap 0.01 \
        --run_ancestry /mnt/nas/Genomics/Genome/FamilyRisk/Datasets/pgsc_HGDP+1kGP_v1.tar.zst > $OUTPUT_DIR/pgscalc.log 2>&1

elif [[ -n "$PGPID" ]]; then
    echo "Use PGPID: $PGPID. "
    $nextflow run pgscatalog/pgsc_calc \
        --input $OUTPUT_DIR/${INPUT_SAMPLE}-pgsc-input.csv \
        --pgp_id $PGPID \
        --target_build $target_build \
        --min_overlap 0.01 \
        --run_ancestry /mnt/nas/Genomics/Genome/FamilyRisk/Datasets/pgsc_HGDP+1kGP_v1.tar.zst > $OUTPUT_DIR/pgscalc.log 2>&1

elif [[ -n "$EFOID" ]]; then
    echo "Use EFOID: $EFOID. "
    $nextflow run pgscatalog/pgsc_calc \
        --input $OUTPUT_DIR/${INPUT_SAMPLE}-pgsc-input.csv \
        --efo_id $EFOID \
        --target_build $target_build \
        --min_overlap 0.01 \
        --run_ancestry /mnt/nas/Genomics/Genome/FamilyRisk/Datasets/pgsc_HGDP+1kGP_v1.tar.zst > $OUTPUT_DIR/pgscalc.log 2>&1

fi

# echo "PRS Pipeline completed at $(date).Please check the results in $OUTPUT_DIR/results."

# ## Add the loop for calculating the other PGS scores directly from the managed vcf files and pfiles. 