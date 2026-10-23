# Load modules
module load SPAdes/4.0.0-foss-2023a-Python-3.11.6

# Paths
INDIR=/nesi/...
OUTDIR=/nesi/...
# SAMPLES=/nesi/.../samples.txt

mkdir -p $OUTDIR

# Get sample ID from list
# SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)
SAMPLE=6109_076_S133

# Input FASTQs
R1=${INDIR}/${SAMPLE}_unmerged_R1.fastq.gz
R2=${INDIR}/${SAMPLE}_unmerged_R2.fastq.gz
MERGED=${INDIR}/${SAMPLE}_merged.fastq.gz

# Output directory per sample
OUT=${OUTDIR}/${SAMPLE}
mkdir -p $OUT

# Run BayesHammer error correction
metaspades.py \
    --only-error-correction \
    -m 340 \
    -k 21,33,55,77 \
    --phred-offset 33 \
    -1 $R1 -2 $R2 \
    --merged $MERGED \
    -o $OUT

