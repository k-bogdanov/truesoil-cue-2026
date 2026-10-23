# Load modules
module purge
module load GTDB-Tk/2.1.1-gimkl-2022a-Python-3.10.5

# File Pathways
IN=/home/~/dRep_output/dereplicated_genomes
OUT=/home/~/gtdbtk_output/

mkdir -p $OUT

gtdbtk classify_wf \
  --genome_dir $IN \
  --out_dir $OUT \
  -x fa \
  --cpus ${SLURM_CPUS_PER_TASK}
