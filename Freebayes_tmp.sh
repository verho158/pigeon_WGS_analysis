#!/bin/bash
#SBATCH --job-name=freebayes_148_samples
#SBATCH --time=9-00:00:00
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=32000
#SBATCH --qos=Std
#SBATCH --output=/lustre/nobackup/WUR/ABGC/verho158/slurm_out_freebayes_%A.txt
#SBATCH --error=/lustre/nobackup/WUR/ABGC/verho158/slurm_err_freebayes_%A.txt

# Set TMPDIR to home directory (writable)
export TMPDIR=$HOME/tmp_freebayes
mkdir -p "$TMPDIR" || { echo "ERROR: TMPDIR not writable"; exit 1; }
echo "Using TMPDIR: $TMPDIR"
df -h "$TMPDIR"

# Load conda environment
source /home/WUR/verho158/miniconda3/etc/profile.d/conda.sh
conda activate pigeon_env

# Input / output directories
bam_list_file="/lustre/nobackup/WUR/ABGC/verho158/all_bams/bam_list.txt"
ref="/lustre/nobackup/WUR/ABGC/verho158/bwa-mem2/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.fna"
scripts_dir="/lustre/nobackup/WUR/ABGC/verho158/scripts"
regions="/lustre/nobackup/WUR/ABGC/verho158/scripts/regions_all_chr.txt"

outdir="/lustre/nobackup/WUR/ABGC/verho158/results/variant_calling"
mkdir -p "$outdir"

vcf_out="${outdir}/pigeon_population_joint_freebayes.vcf.gz"

echo "[$(date)] Starting joint Freebayes calling on 148 samples" >&2
echo "Using BAM list: $bam_list_file" >&2

# Run Freebayes joint calling
${scripts_dir}/freebayes-parallel.sh \
    ${regions} 12 \
    -f "$ref" \
    --use-best-n-alleles 2 \
    --min-base-quality 10 \
    --min-alternate-fraction 0.2 \
    --min-alternate-count 2 \
    --haplotype-length 0 \
    --ploidy 2 \
    -L "$bam_list_file" \
| vcffilter -f "QUAL > 20 & DP > 3" \
| bgzip -c > "$vcf_out"

# Index the VCF
tabix -p vcf "$vcf_out"

echo "[$(date)] Finished joint Freebayes calling" >&2

conda deactivate
