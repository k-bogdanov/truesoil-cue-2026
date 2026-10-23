
# Modules
module load Java
module load BBMap
module load Trimmomatic/0.39-Java-1.8.0_144

# Paths
HG38=/nesi/.../hg38.fa
INDIR=/nesi/...
OUTDIR=/nesi/...
SAMPLES=/nesi/.../samples.txt

# Get sample name
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)

# Human contaminant removal (BBMap)
bbmap.sh ref=$HG38 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
    qtrim=rl trimq=10 untrim -Xmx50g \
    in1=${OUTDIR}/${SAMPLE}_R1_clean.fastq.gz \
    in2=${OUTDIR}/${SAMPLE}_R2_clean.fastq.gz \
    outu1=${OUTDIR}/${SAMPLE}_R1_hostRemoved.fastq.gz \
    outu2=${OUTDIR}/${SAMPLE}_R2_hostRemoved.fastq.gz \
    outm1=${OUTDIR}/${SAMPLE}_R1_human.fastq.gz \
    outm2=${OUTDIR}/${SAMPLE}_R2_human.fastq.gz

# BBMerge
bbmerge-auto.sh \
    in1=${OUTDIR}/${SAMPLE}_R1_hostRemoved.fastq.gz \
    in2=${OUTDIR}/${SAMPLE}_R2_hostRemoved.fastq.gz \
    out=${OUTDIR}/${SAMPLE}_merged.fastq.gz \
    outu1=${OUTDIR}/${SAMPLE}_unmerged_R1.fastq.gz \
    outu2=${OUTDIR}/${SAMPLE}_unmerged_R2.fastq.gz \
    t=16 prealloc rem k=5 extend2=50 ecct -Xmx48g

