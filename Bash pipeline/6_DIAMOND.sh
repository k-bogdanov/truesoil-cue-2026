# Load modules
module load DIAMOND/2.1.14-GCC-12.3.0

# File pathways
SAMPLES=/nesi/.../samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES")
DB_FOLDER=/nesi/.../Tools/greening/greening_db_dmnd
INDIR=/nesi/.../${SAMPLE}/prodigal
OUTDIR=/nesi/.../${SAMPLE}/greening_hits
mkdir -p "$OUTDIR"

# Diamond
for DB in "$DB_FOLDER"/*.dmnd; do
    NAME=$(basename "$DB" .dmnd)
    diamond blastp \
      --query "$INDIR/proteins.faa" \
      --db "$DB" \
      --out "$OUTDIR/${NAME}_hits.tsv" \
      --outfmt 6 qseqid sseqid pident length evalue bitscore \
      --evalue 1e-5 --max-target-seqs 1 --threads 4
done

