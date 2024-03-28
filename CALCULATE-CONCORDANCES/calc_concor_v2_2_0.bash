#!/bin/bash
#SBATCH --job-name=CONCOR22
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bj2349@nyulangone.org
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --output=calc_concor_v2_2_0_%j.log
#SBATCH -p cpu_short

### [[Welcome]] ###

echo "-------------------------------------------------------------"
echo "Welcome to the calc_concordance tool!"
echo "Functionality: Calculate genotype concordance and non-reference genotype concordance"
echo "Version: V2.2.0"
echo "Copyright: Baiqi Jiang & Chakravarti Lab @2024"
echo "For more information, please refer to readme.txt"
echo "-------------------------------------------------------------"

### [[Initialization]] ###

echo "Initializing parameters..."

# Parameters
compare_sample="HG00099alt2.chr21.vcf"
true_sample="HG00099ori.chr21.vcf"
target_bin_info="AF"  # or TAF, or any other allele frequency info field
reference_file="/gpfs/scratch/bj2349/Research_W5_Mar17_23/reference/GRCh38_latest_genomic.fna"
delete_intermediate_files=true
add_reference_file=false

# Validate reference file parameter
if $add_reference_file && [ -z "$reference_file" ]; then
    echo "When add_reference_file is set to true, you must specify reference_file."
    exit 1
fi

# Conditional reference parameter
reference_param=""
if $add_reference_file; then
    reference_param="-R $reference_file"
fi

# Remove file extensions for directory names
compare_sample_name="${compare_sample%%.*}"
true_sample_name="${true_sample%%.*}"

# Directories
WORK_DIR="/gpfs/scratch/bj2349/Research_W6_Mar24_30"
COMPARE_DIR="${WORK_DIR}/${compare_sample_name}"
TRUE_DIR="${WORK_DIR}/${true_sample_name}"
GATK_CONCORDANCE_DIR="${WORK_DIR}/${compare_sample_name}_vs_${true_sample_name}_concordance_GATK"
CONCORDANCE_DIR="${WORK_DIR}/${compare_sample_name}_vs_${true_sample_name}_concordance"

echo "Parameters initialized."

### [[Module Loading]] ###

echo "Loading modules..."
module load gatk
module load bcftools
echo "Modules are loaded."

### [[Pre-set Check and Indexing]] ###

# Check and index if needed
echo "Checking and indexing VCF files..."
if [ ! -f "${WORK_DIR}/${compare_sample}.csi" ]; then
    echo "Indexing ${compare_sample}..."
    bcftools index "${WORK_DIR}/${compare_sample}"
fi
if [ ! -f "${WORK_DIR}/${true_sample}.csi" ]; then
    echo "Indexing ${true_sample}..."
    bcftools index "${WORK_DIR}/${true_sample}"
fi

### [[Directories Creation]] ###

echo "Creating directories..."
mkdir -p "${COMPARE_DIR}"
mkdir -p "${TRUE_DIR}"
mkdir -p "${GATK_CONCORDANCE_DIR}"
mkdir -p "${CONCORDANCE_DIR}"
echo "Directories created."

### [[Intervals Handling]] ###

# Define intervals based on target_bin_info
declare -a intervals=("0-0.01" "0.01-0.02" "0.02-0.05" "0.05-0.1" "0.1-0.2" "0.2-0.5" "0.5-1")

# Process intervals
for interval in "${intervals[@]}"; do
    echo "Processing interval: ${interval}"
    IFS='-' read -r start end <<< "$interval"

    # Filtering and indexing
    echo "Filtering based on ${target_bin_info} for ${interval}..."
    bcftools view -i "INFO/${target_bin_info}>=${start} && INFO/${target_bin_info}<=${end}" "${WORK_DIR}/${compare_sample}" -Oz -o "${COMPARE_DIR}/${compare_sample_name}_${target_bin_info}_${interval}.vcf.gz"
    bcftools view -i "INFO/${target_bin_info}>=${start} && INFO/${target_bin_info}<=${end}" "${WORK_DIR}/${true_sample}" -Oz -o "${TRUE_DIR}/${true_sample_name}_${target_bin_info}_${interval}.vcf.gz"
    bcftools index "${COMPARE_DIR}/${compare_sample_name}_${target_bin_info}_${interval}.vcf.gz"
    bcftools index "${TRUE_DIR}/${true_sample_name}_${target_bin_info}_${interval}.vcf.gz"

    # GenotypeConcordance
    echo "Calculating GenotypeConcordance for ${interval}..."
    gatk GenotypeConcordance \
        $reference_param \
        --CALL_VCF "${COMPARE_DIR}/${compare_sample_name}_${target_bin_info}_${interval}.vcf.gz" \
        --CALL_SAMPLE "${compare_sample_name}" \
        --TRUTH_VCF "${TRUE_DIR}/${true_sample_name}_${target_bin_info}_${interval}.vcf.gz" \
        --TRUTH_SAMPLE "${true_sample_name}" \
        --OUTPUT_ALL_ROWS \
        -O "${GATK_CONCORDANCE_DIR}/gc_results_${compare_sample_name}_vs_${true_sample_name}_${target_bin_info}_${interval}"
done

# GATK checker

if [ "$(ls -A $GATK_CONCORDANCE_DIR)" ]; then
    echo "GATK has successfully run across all user-specified INFO intervals, with all files stored in $GATK_CONCORDANCE_DIR."
else
    echo "FAILED: The $GATK_CONCORDANCE_DIR directory is empty, indicating a potential failure in the GATK run. For more information, please refer back to the section on running GATK."
fi

### [[Calculate Concordance and Save to .tsv]] ###

echo "Calculating concordance metrics..."

# Create output file for concordance metrics
FINAL_TSV="${CONCORDANCE_DIR}/gc_results_${compare_sample_name}_vs_${true_sample_name}_all_intervals.tsv"
echo -e "VAR_TYPE\t${target_bin_info}_INTERVAL\tG_CONCORDANCE\tNON_REF_G_CONCORDANCE" > "$FINAL_TSV"

echo "Calculating concordance metrics for each interval and saving to a single TSV file..."

# Loop through the intervals and calculate concordance
for interval in "${intervals[@]}"; do
    METRICS_FILE="${GATK_CONCORDANCE_DIR}/gc_results_${compare_sample_name}_vs_${true_sample_name}_${target_bin_info}_${interval}.genotype_concordance_contingency_metrics"
    
    # Extract TP, TN, FP, FN counts for SNP and INDEL from the contingency metrics file
    while IFS=$'\t' read -r -a line; do
        if [[ "${line[0]}" == "SNP" || "${line[0]}" == "INDEL" ]]; then
            TP="${line[3]}"
            TN="${line[4]}"
            FP="${line[5]}"
            FN="${line[6]}"
            
            # Calculate GENOTYPE_CONCORDANCE and NON_REF_GENOTYPE_CONCORDANCE
            G_CONCORDANCE=$(echo "scale=4; ($TP + $TN) / ($TP + $TN + $FP + $FN)" | bc)
            NON_REF_G_CONCORDANCE=$(echo "scale=4; $TP / ($TP + $FP + $FN)" | bc)
            
            # Append to the final .tsv file
            echo -e "${line[0]}\t${interval}\t$G_CONCORDANCE\t$NON_REF_G_CONCORDANCE" >> "$FINAL_TSV"
        fi
    done < <(tail -n +8 "$METRICS_FILE") # Skip headers
    
    echo "Concordance metrics for interval ${interval} appended to $FINAL_TSV"
done

# Final .tsv checker

if [ "$(ls -A $CONCORDANCE_DIR)" ] && [ -s "$FINAL_TSV" ]; then
    echo "All intervals processed and concordance metrics saved to $FINAL_TSV"
elif [ -f "$FINAL_TSV" ] && ! [ -s "$FINAL_TSV" ]; then
    echo "WARNING: $FINAL_TSV was successfully generated, but the file is empty. Please check the files with the suffix 'contingency matrices' in the $GATK_CONCORDANCE_DIR, as well as the version of GATK."
else
    echo "FAILED: Concordance calculation failed, no files were generated. Please check the preliminary steps."
fi

### [[Cleanup]] ###

if [ "$delete_intermediate_files" = "true" ]; then
    echo "Deleting intermediate directories..."
    rm -r "${COMPARE_DIR}"
    rm -r "${TRUE_DIR}"
    echo "Intermediate directories deleted."
fi
