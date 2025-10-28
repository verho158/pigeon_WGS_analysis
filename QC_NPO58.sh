#!/bin/bash
#SBATCH --job-name=QC_NPO58
#SBATCH --time=1-00:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16000
#SBATCH --qos=Std
#SBATCH --output=/lustre/nobackup/WUR/ABGC/verho158/QC/slurm_output_QC_NPO58_%j.txt
#SBATCH --error=/lustre/nobackup/WUR/ABGC/verho158/QC/slurm_error_QC_NPO58_%j.txt

# Variables
sample="NPO_58"
reads_dir="/lustre/nobackup/WUR/ABGC/verho158/WGS_Data/${sample}"
out_dir="/lustre/nobackup/WUR/ABGC/verho158/QC"

# Create output directories if they don't exist
mkdir -p $out_dir
filtered_dir="${out_dir}/filtered"
mkdir -p $filtered_dir

# Load Conda environment
source /home/WUR/verho158/miniconda3/etc/profile.d/conda.sh
conda activate pigeon_env

# Step 1: Quality check (FastQC)
echo "Running FastQC for ${sample}..."
fastqc -t 4 -o $out_dir \
  ${reads_dir}/${sample}_EKDN250028352-1A_22WWJCLT4_L4_1.fq.gz \
  ${reads_dir}/${sample}_EKDN250028352-1A_22WWJCLT4_L4_2.fq.gz

# Step 2: Filtering and trimming (fastp)
echo "Running fastp filtering for ${sample}..."
fastp \
  -i ${reads_dir}/${sample}_EKDN250028352-1A_22WWJCLT4_L4_1.fq.gz \
  -I ${reads_dir}/${sample}_EKDN250028352-1A_22WWJCLT4_L4_2.fq.gz \
  -o ${filtered_dir}/${sample}_filtered_1.fq.gz \
  -O ${filtered_dir}/${sample}_filtered_2.fq.gz \
  -q 30 \
  -l 36 \
  -w 4 \
  -D \
  -h ${filtered_dir}/${sample}_fastp.html \
  -j ${filtered_dir}/${sample}_fastp.json

# Step 3: Summarize with MultiQC
echo "Summarizing QC results with MultiQC..."
multiqc $out_dir -o $out_dir

# Clean up
conda deactivate
echo "QC and filtering completed successfully for ${sample}!"
