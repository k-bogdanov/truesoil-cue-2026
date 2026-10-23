# Load modules
module purge
module load metaWRAP/1.2.1-Miniconda3

# File Pathways
SAMPLES=/nesi/~/samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)

# Input dirs
AS=/home/~/${SAMPLE}/4_spades/spades
READS=/home/~/${SAMPLE}/3_errors_corrected/corrected

OUT=/nesi/~/${SAMPLE}/metaWRAP_bins
mkdir -p $OUT

metawrap binning -o $OUT \
                 -t $SLURM_CPUS_PER_TASK \
                 -a ${AS}/contigs.fasta \
                 --metabat2 --maxbin2 --concoct \
                 ${READS}/${SAMPLE}_unmerged_1.fastq \
                 ${READS}/${SAMPLE}_unmerged_2.fastq