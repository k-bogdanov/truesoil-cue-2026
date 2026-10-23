# Load modules
module load SPAdes/4.0.0-foss-2023a-Python-3.11.6

# Paths
SAMPLES=/nesi/.../samples.txt

# Get sample name
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)

INDIR=/nesi/...
OUTDIR=/nesi/...

# Input and output paths
IN=${INDIR}/${SAMPLE}/corrected
OUT=${OUTDIR}/${SAMPLE}
mkdir -p $OUT

# Run assembly
metaspades.py \
    --only-assembler \
    -m 340 \
    -k 21,33,55,77 \
    --threads 16 \
    --phred-offset 33 \
    -1 $IN/${SAMPLE}_unmerged_R1.fastq00.0_0.cor.fastq.gz \
    -2 $IN/${SAMPLE}_unmerged_R2.fastq00.0_0.cor.fastq.gz \
    --merged $IN/${SAMPLE}_merged.fastq00.0_1.cor.fastq.gz \
    -s $IN/${SAMPLE}_unmerged_R_unpaired00.0_0.cor.fastq.gz \
    -o $OUT/spades
