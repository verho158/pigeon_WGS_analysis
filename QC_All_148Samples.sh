#!/bin/bash
#SBATCH --job-name=QC_NPO
#SBATCH --time=1-00:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16000
#SBATCH --qos=Std
#SBATCH --output=/lustre/nobackup/WUR/ABGC/verho158/QC/slurm_output_QC_%A_%a.txt
#SBATCH --error=/lustre/nobackup/WUR/ABGC/verho158/QC/slurm_error_QC_%A_%a.txt
#SBATCH --array=1-148 # Default full array (override in sbatch command for testing)

# Load Conda environment
source /home/WUR/verho158/miniconda3/etc/profile.d/conda.sh
conda activate pigeon_env

# Base directories
base_dir="/lustre/nobackup/WUR/ABGC/verho158/WGS_Data"
out_dir="/lustre/nobackup/WUR/ABGC/verho158/QC"

# Determine sample name
sample=$(printf "NPO_%d" $SLURM_ARRAY_TASK_ID)
reads_dir="${base_dir}/${sample}"

# Create sample-specific directories
sample_out_dir="${out_dir}/${sample}"
filtered_dir="${sample_out_dir}/filtered"
mkdir -p $sample_out_dir $filtered_dir

# Find FASTQ files
fastq_files=(${reads_dir}/*.fq.gz)
if [ ${#fastq_files[@]} -eq 0 ]; then
    echo "No FASTQ files found for ${sample}, skipping..."
    exit 1
fi

# Step 1: Run FastQC
echo "Running FastQC for ${sample}..."
fastqc -t 4 -o $sample_out_dir "${fastq_files[@]}"

# Step 2: Run fastp filtering
echo "Running fastp for ${sample}..."
if [ ${#fastq_files[@]} -eq 2 ]; then
    fastp \
      -i ${fastq_files[0]} \
      -I ${fastq_files[1]} \
      -o ${filtered_dir}/${sample}_filtered_1.fq.gz \
      -O ${filtered_dir}/${sample}_filtered_2.fq.gz \
      -q 30 -l 36 -w 4 -D \
      -h ${filtered_dir}/${sample}_fastp.html \
      -j ${filtered_dir}/${sample}_fastp.json
else
    # Handle multiple paired FASTQs
    for ((i=0; i<${#fastq_files[@]}; i+=2)); do
        base_name="${sample}_filtered_$((i/2+1))"
        fastp \
          -i ${fastq_files[i]} \
          -I ${fastq_files[i+1]} \
          -o ${filtered_dir}/${base_name}_1.fq.gz \
          -O ${filtered_dir}/${base_name}_2.fq.gz \
          -q 30 -l 36 -w 4 -D \
          -h ${filtered_dir}/${base_name}_fastp.html \
          -j ${filtered_dir}/${base_name}_fastp.json
    done
fi

# Deactivate Conda
conda deactivate
echo "QC and filtering completed successfully for ${sample}!"
