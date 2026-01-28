#!/bin/bash
#SBATCH --job-name=bwamem2_index
#SBATCH -c 4
#SBATCH --mem=16000
#SBATCH --time=02:00:00
#SBATCH --output=/lustre/nobackup/WUR/ABGC/verho158/bwamem2_index_%j.txt
#SBATCH --error=/lustre/nobackup/WUR/ABGC/verho158/bwamem2_index_%j.err

# Load Conda environment
source /home/WUR/verho158/miniconda3/etc/profile.d/conda.sh
conda activate pigeon_env

# Reference genome path
ref="/lustre/nobackup/WUR/ABGC/verho158/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.fna"

# Run BWA-MEM2 indexing
echo "Indexing reference genome with BWA-MEM2..."
bwa-mem2 index $ref

echo "Indexing complete."
conda deactivate
