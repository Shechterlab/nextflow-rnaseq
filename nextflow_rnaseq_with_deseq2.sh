#!/bin/bash
#This script uses Nextflow to process RNA-seq data with nextflow then run DESEQ2 on the output
#This script requires a sample metadata file called "samplesheet.csv"
#e.g. Run as follows: "sbatch nextflow_rnaseq_with_deseq2.sh"

#By Maxim Maron on 12-09-22

#SBATCH -p unlimited       #partition/queue name
#SBATCH --job-name=nextflow_rnaseq    # Job name
#SBATCH --mail-type=END,FAIL    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maxim.maron@einsteinmed.edu # Where to send mail
#SBATCH --ntasks=1      # Run on a single CPU
#SBATCH --cpus-per-task=40      # Number of CPU cores per task
#SBATCH --mem=80gb      # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=nextflow_RNAseq.log      # Standard output and error log

#Load conda envrionment and modules
source  /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate
conda activate zlib_nextflow_rnaseq
module load singularity
module load nextflow

#Export source variables to envrionment
.  /gs/gsfs0/users/shechter-lab/data/NGS/stds/source_files/20221210_nextflow_source_file.sh

nextflow run nf-core/rnaseq --input samplesheet.csv \
-profile singularity \
--outdir ./nextflow_results \
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

mkdir ./nextflow_results/DESEQ2

#Need to make DESEQ2 metadata sheet
awk -f $reformat_metadata_file samplesheet.csv > ./nextflow_results/DESEQ2/samplesheet_DESeq2.csv

conda activate R_nextflow_rnaseq

#Execute R script with reference files in positional arguments from source variables
Rscript $DESEQ2 $tx2gene $gtf_gene_annotation_table

conda deactivate

#Test if run is completed successfully by checking log file and if so delete the work directory
 if grep -q "Pipeline completed successfully" nextflow_RNAseq.log; then   
 rm -r ./work; 
 fi

