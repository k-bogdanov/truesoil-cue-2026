# Load modules
module load prodigal/2.6.3-GCCcore-7.4.0

# File pathways
SAMPLES=/nesi/.../samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)
INDIR=/nesi/.../${SAMPLE}/spades
OUTDIR=/nesi/.../${SAMPLE}/prodigal
mkdir -p $OUTDIR

# Prodigal
cd $OUTDIR

prodigal \
  -i $INDIR/contigs.fasta \
  -a proteins.faa \
  -d genes.fna \
  -o prodigal.txt \
  -p meta -m -q

