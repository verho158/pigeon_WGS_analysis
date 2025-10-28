#!/bin/bash
#SBATCH --job-name=mapping_NPO58
#SBATCH --time=2-00:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16000
#SBATCH --qos=Std
#SBATCH --output=/lustre/nobackup/WUR/ABGC/verho158/mapping/slurm_output_mapping_NPO58_%j.txt
#SBATCH --error=/lustre/nobackup/WUR/ABGC/verho158/mapping/slurm_error_mapping_NPO58_%j.txt

# Activate Conda environment
source /home/WUR/verho158/miniconda3/etc/profile.d/conda.sh
conda activate pigeon_env

# Variables
sample="NPO_58"
ref="/lustre/nobackup/WUR/ABGC/verho158/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.fna"
reads_dir="/lustre/nobackup/WUR/ABGC/verho158/QC/filtered"
out_dir="/lustre/nobackup/WUR/ABGC/verho158/mapping/${sample}"

mkdir -p $out_dir

# Step 1: Mapping with BWA-MEM
echo "Mapping ${sample} reads to reference genome..."
bwa mem -t 4 $ref \
  ${reads_dir}/${sample}_filtered_1.fq.gz \
  ${reads_dir}/${sample}_filtered_2.fq.gz \
  > ${out_dir}/${sample}.sam

# Step 2: Convert, sort, and index BAM file
echo "Converting SAM to sorted BAM..."
samtools view -@ 4 -bS ${out_dir}/${sample}.sam | samtools sort -@ 4 -o ${out_dir}/${sample}_sorted.bam
samtools index ${out_dir}/${sample}_sorted.bam

# Step 3: Generate mapping statistics
echo "Generating alignment stats..."
samtools flagstat ${out_dir}/${sample}_sorted.bam > ${out_dir}/${sample}_flagstat.txt

# Step 4: Calculate average coverage
echo "Calculating coverage..."
samtools depth -a ${out_dir}/${sample}_sorted.bam | \
  awk '{sum+=$3; count++} END {if (count>0) print "Average coverage:", sum/count; else print "No data."}' \
  > ${out_dir}/${sample}_coverage.txt

# Deactivate Conda
conda deactivate

echo "Mapping and coverage analysis completed successfully for ${sample}!"
