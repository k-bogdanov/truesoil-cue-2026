# Header example:
#!/bin/bash
#SBATCH --account=#account
#SBATCH --job-name=trimm
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-5              
#SBATCH --output /nesi/.../logs/%x_%A_%a.out
#SBATCH --error /nesi/.../logs/%x_%A_%a.err

# Modules
module load Java
module load BBMap
module load Trimmomatic/0.39-Java-1.8.0_144

# Paths
INDIR=/nesi/...
OUTDIR=/nesi/...
ADAPTERS=/nesi/.../bbmap/NexteraPE-PE.fa
SAMPLES=/nesi/.../samples.txt

# Get sample name
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)

R1=${INDIR}/${SAMPLE}_R1_001.fastq.gz
R2=${INDIR}/${SAMPLE}_R2_001.fastq.gz

# Clumpify
clumpify.sh in1=$R1 in2=$R2 \
    out1=${OUTDIR}/${SAMPLE}_R1_clust.fastq.gz \
    out2=${OUTDIR}/${SAMPLE}_R2_clust.fastq.gz

# Trimmomatic
trimmomatic PE -threads 16 -phred33 \
    ${OUTDIR}/${SAMPLE}_R1_clust.fastq.gz ${OUTDIR}/${SAMPLE}_R2_clust.fastq.gz \
    ${OUTDIR}/${SAMPLE}_R1_clean.fastq.gz ${OUTDIR}/${SAMPLE}_R1_unpaired.fastq.gz \
    ${OUTDIR}/${SAMPLE}_R2_clean.fastq.gz ${OUTDIR}/${SAMPLE}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:1:25:7 SLIDINGWINDOW:4:30 MINLEN:80

