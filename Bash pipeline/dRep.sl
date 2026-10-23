# Load dRep
module load drep/3.4.2-gimkl-2022a-Python-3.10.5

# Pathways

IN=/home/~/refined_bins
OUT=/home/~/dRep_output
mkdir -p $OUT

dRep dereplicate $OUT -g ${IN}/*.fa --ignoreGenomeQuality --processors 16 --skip_plots --S_algorithm fastANI