#!/bin/bash
#test nextflow rnaseq pipeline

#SBATCH -p unlimited       #partition/queue name
#SBATCH --job-name=nextflow_rna    # Job name
#SBATCH --mail-type=END,FAIL    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maxim.maron@einsteinmed.edu # Where to send mail
#SBATCH --ntasks=1      # Run on a single CPU
#SBATCH --cpus-per-task=40      # Number of CPU cores per task
#SBATCH --mem=80gb      # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=nextflow_RNAseq_test.log      # Standard output and error log

source  /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate
conda activate zlib_nextflow_rnaseq
module load singularity
module load nextflow

nextflow run nf-core/rnaseq -profile test,singularity --outdir nextflow_test

module unload nextflow
moedule unload singularity
conda deactivate
