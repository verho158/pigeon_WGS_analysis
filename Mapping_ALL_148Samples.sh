#!/bin/bash
#SBATCH --job-name=map_NPO
#SBATCH --time=2-00:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16000
#SBATCH --qos=Std
#SBATCH --output=/lustre/nobackup/WUR/ABGC/verho158/mapping/slurm_out_map_%A_%a.txt
#SBATCH --error=/lustre/nobackup/WUR/ABGC/verho158/mapping/slurm_err_map_%A_%a.txt
#SBATCH --array=1-148% # Default full array (override in sbatch command for testing)

# ---- load env ----
source /home/WUR/verho158/miniconda3/etc/profile.d/conda.sh
conda activate pigeon_env

# ---- vars ----
sample=$(printf "NPO_%d" $SLURM_ARRAY_TASK_ID)
reads_dir="/lustre/nobackup/WUR/ABGC/verho158/QC/${sample}/filtered"
out_dir="/lustre/nobackup/WUR/ABGC/verho158/mapping/${sample}"
ref="/lustre/nobackup/WUR/ABGC/verho158/bwa-mem2/GCF_036013475.1_bColLiv1.pat.W.v2_genomic.fna"

mkdir -p "$out_dir"

# ---- find read pairs (R1/R2 pattern) ----
R1=( ${reads_dir}/*_1.fq.gz )
R2=( ${reads_dir}/*_2.fq.gz )

if [ ${#R1[@]} -eq 0 ] || [ ${#R1[@]} -ne ${#R2[@]} ]; then
  echo "ERROR: no paired FASTQs or unmatched pairs for ${sample} in ${reads_dir}" >&2
  exit 1
fi

# ---- map each pair to temporary BAM part files ----
parts=()
for ((i=0; i<${#R1[@]}; i++)); do
  r1=${R1[i]}
  r2=${R2[i]}
  part="${out_dir}/${sample}_part_$((i+1)).bam"
  echo "[$(date)] Mapping pair $((i+1)) for ${sample}: $r1 + $r2" >&2
  bwa-mem2 mem -t ${SLURM_CPUS_ON_NODE:-4} "$ref" "$r1" "$r2" \
    | samtools view -@ ${SLURM_CPUS_ON_NODE:-4} -b -o "$part" -
  parts+=( "$part" )
done

# ---- if multiple parts, merge, else use single part ----
if [ ${#parts[@]} -gt 1 ]; then
  merged="${out_dir}/${sample}_merged.bam"
  echo "[$(date)] Merging ${#parts[@]} parts into $merged" >&2
  samtools merge -@ ${SLURM_CPUS_ON_NODE:-4} -f "$merged" "${parts[@]}"
  rm -f "${parts[@]}"
  input_bam="$merged"
else
  input_bam="${parts[0]}"
fi

# ---- sort (coordinate) ----
sorted="${out_dir}/${sample}_sorted.bam"
echo "[$(date)] Sorting to $sorted" >&2
samtools sort -@ ${SLURM_CPUS_ON_NODE:-4} -o "$sorted" "$input_bam"
rm -f "$input_bam"

# ---- add read groups (using Picard) ----
RGID="${sample}"
RGLB="lib1"
RGPL="ILLUMINA"
RGPU="unit1"
RGSM="${sample}"

rg_bam="${out_dir}/${sample}_RG.bam"
echo "[$(date)] Adding read groups" >&2
picard AddOrReplaceReadGroups \
  I="$sorted" O="$rg_bam" \
  RGID="$RGID" RGLB="$RGLB" RGPL="$RGPL" RGPU="$RGPU" RGSM="$RGSM"
rm -f "$sorted"

# ---- mark duplicates (Picard) and create index ----
dedup="${out_dir}/${sample}_dedup.bam"
dupmetrics="${out_dir}/${sample}_dupmetrics.txt"
echo "[$(date)] Marking duplicates" >&2
picard MarkDuplicates I="$rg_bam" O="$dedup" M="$dupmetrics" CREATE_INDEX=true
rm -f "$rg_bam"

# ---- flagstat and coverage summary ----
samtools flagstat "$dedup" > "${out_dir}/${sample}_flagstat.txt"

samtools depth -a "$dedup" \
  | awk '{sum+=$3; if($3>=1) c1++; if($3>=5) c5++; if($3>=10) c10++; if($3>=20) c20++; count++}
         END {if(count>0) {print "avg:",sum/count; printf ">=1x: %.2f%%\n>=5x: %.2f%%\n>=10x: %.2f%%\n>=20x: %.2f%%\n",
               100*c1/count, 100*c5/count, 100*c10/count, 100*c20/count } else print "No coverage"}' \
  > "${out_dir}/${sample}_depth_summary.txt"

# ---- done ----
conda deactivate
echo "[$(date)] Sample ${sample} mapping DONE" >&2
