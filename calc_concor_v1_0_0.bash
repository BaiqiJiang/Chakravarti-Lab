#!/bin/bash
#SBATCH --job-name=CONCOR
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bj2349@nyulangone.org
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --output=calc_concor_v1_0_0_%j.log
#SBATCH -p cpu_short

### [[Initialization]] ###

echo "Initializing parameters..."

# Parameters
compare_sample="HG00099.chr21.vcf.gz"
true_sample="HG00117.chr21.vcf.gz"
target_bin_info="AF"  # or TAF, or any other allele frequency info field
reference_file="/gpfs/scratch/bj2349/Research_W4_Mar10_16/reference/GRCh38_latest_genomic.fna"
delete_intermediate_files=true

# Remove file extensions for directory names
compare_sample_name="${compare_sample%%.*}"
true_sample_name="${true_sample%%.*}"

# Directories
WORK_DIR="/gpfs/scratch/bj2349/Research_W4_Mar10_16"
COMPARE_DIR="${WORK_DIR}/${compare_sample_name}"
TRUE_DIR="${WORK_DIR}/${true_sample_name}"
CONCORDANCE_DIR="${WORK_DIR}/${compare_sample_name}_vs_${true_sample_name}_concordance_ALL_ROWS"

echo "Parameters initialized."

### [[Module Loading]] ###

echo "Loading modules..."
module load gatk
module load bcftools
echo "Modules loaded: GATK, bcftools"

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
        -R "${reference_file}" \
        --CALL_VCF "${COMPARE_DIR}/${compare_sample_name}_${target_bin_info}_${interval}.vcf.gz" \
        --CALL_SAMPLE "${compare_sample_name}" \
        --TRUTH_VCF "${TRUE_DIR}/${true_sample_name}_${target_bin_info}_${interval}.vcf.gz" \
        --TRUTH_SAMPLE "${true_sample_name}" \
	--OUTPUT_ALL_ROWS \
        -O "${CONCORDANCE_DIR}/gc_results_${compare_sample_name}_vs_${true_sample_name}_${interval}"
done

### [[Cleanup]] ###

if [ "$delete_intermediate_files" = "true" ]; then
    echo "Deleting intermediate directories..."
    rm -r "${COMPARE_DIR}"
    rm -r "${TRUE_DIR}"
    echo "Intermediate directories deleted."
fi

echo "All tasks completed."
