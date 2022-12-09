#!/bin/bash
#This script uses Nextflow to process RNA-seq data with nextflow then run DESEQ2 on the output
#This script requires a sample metadata file called "nextflow_sample_metadata.csv" that is located in the premade output directory
#Only one argument is necessary: full path to desired output directory with no trailing forward slash
#e.g. Run as follows: "sbatch nextflow_rnaseq_with_deseq2.sh /path/to/output/dir"

#By Maxim Maron on 12-09-22

#SBATCH -p unlimited       #partition/queue name
#SBATCH --job-name=nextflow_rnaseq    # Job name
#SBATCH --mail-type=END,FAIL    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maxim.maron@einsteinmed.edu # Where to send mail
#SBATCH --ntasks=1      # Run on a single CPU
#SBATCH --cpus-per-task=40      # Number of CPU cores per task
#SBATCH --mem=80gb      # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=nextflow_RNA_test.log      # Standard output and error log

#Load conda envrionment and modules
source  /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate
conda activate zlib_nextflow_rnaseq
module load singularity
module load nextflow

#Export source variables to envrionment
.  /gs/gsfs0/users/shechter-lab/data/NGS/stds/source_files/20221207_nextflow_source_file.sh

nextflow run nf-core/rnaseq --input $1/nextflow_sample_metadata.csv \
-profile singularity \
--outdir $1 \
--aligner star_salmon \
--fasta $nextflow_rna_seq_fasta \
--gtf $nextflow_rna_seq_gtf \
--star_index $nextflow_rna_seq_star_index \
--rsem_index $nextflow_rna_seq_rsem_index \
--gene_bed $nextflow_rna_seq_gene_bed \
--transcript_fasta $nextflow_rna_seq_transcript_fasta

module unload nextflow
module unload singularity
conda deactivate

mkdir $1/DESEQ2

#Need to make DESEQ2 metadata sheet
awk -f $reformat_metadata_file $1/nextflow_sample_metadata.csv > $1/DESEQ2/nextflow_sample_metadata.csv.DESEQ2.metadata.csv

conda activate R_nextflow_rnaseq

Rscript $DESEQ2 $1

conda deactivate
