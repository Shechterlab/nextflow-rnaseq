#!/bin/bash
#This script uses Nextflow to process RNA-seq data with nextflow then run DESEQ2 on the output
#This script requires a sample metadata file called "samplesheet.csv"
#e.g. Run as follows: "sbatch nextflow_rnaseq_with_deseq2.sh"

#By Maxim Maron on 12-09-22

#SBATCH -p unlimited       #partition/queue name
#SBATCH --job-name=nextflow_rnaseq    # Job name
#SBATCH --mail-type=NONE    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maxim.maron@einsteinmed.edu # Where to send mail
#SBATCH --ntasks=1      # Run on a single CPU
#SBATCH --cpus-per-task=40      # Number of CPU cores per task
#SBATCH --mem=80gb      # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=nextflow_RNAseq.log      # Standard output and error log


###RUNNING NEXTFLOW###

#Load conda envrionment and modules
source  /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate
conda activate zlib_nextflow_rnaseq
module load singularity
module load nextflow

#Export source variables to envrionment
.  /gs/gsfs0/users/shechter-lab/data/NGS/stds/source_files/20221220_nextflow_source_file.sh

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

###END NEXTFLOW###

##RUNNING RMATS###
#activate the rmats conda envrionment
conda activate rmats_nextflow

#make a directory to store rmats files
mkdir rmats/

#create an empty array called length that will store read lengths used to determine median read length for rmats
length=()

#go through all of the bam files and determine the read length, then append this value to the length array
#encountered error: "samtools: error while loading shared libraries: libncurses.so.5: cannot open shared object file: No such file or directory" that was fixed by going into the envrionment libraries and renaming libncurses.so.6 to libncurses.so.5

for file in $(ls nextflow_results/star_salmon/*.markdup.sorted.bam); do
	length+=( $(samtools view $file | head -n 1000 | gawk '{print length($10)}' | sort | uniq -c | perl -ane '$_ =~ s/^[ ]+//g;print $_' | sort -k 1nr,1nr | head -1 | cut -f2 -d " ") )
done

#determine the median read length and store this as a variable to be used with rmats
read_length=$(printf "%s\n" "${length[@]}" | datamash median 1)

#Make samplesheet that points to location of bam files for rmats
variables=$(awk -F',' 'FNR>1 && !a[$1]++{print $1}' samplesheet.csv)
for variable in $variables; do
	printf `pwd`/nextflow_results/star_salmon/${variable}.markdup.sorted.bam, >> rmats/${variable%%_REP*}_rmats_tmp.csv
	sed 's/\(.*\),/\1 /' rmats/${variable%%_REP*}_rmats_tmp.csv > rmats/${variable%%_REP*}_rmats.csv
done
rm rmats/*_rmats_tmp.csv


#build an array of file locations to analyze with rmats and remove the CONTROL
files_to_analyze=()
files=$(ls ./rmats/*.csv)
for file in $files; do
 files_to_analyze+=$file
done
delete=./rmats/CONTROL_rmats.csv
files_to_analyze=( "${files_to_analyze[@]/$delete}")

#run rmats using the array of variables
for file in $files_to_analyze; do
#only take the name after the last forward slash
 name=$(echo $file | sed 's:.*/::')
 rmats.py --b1 ./rmats/CONTROL_rmats.csv --b2 $file --gtf $rmats_gtf -t paired --readLength $read_length --nthread 40 --od ./rmats/CONTROL_vs_${name%_rmats.csv} --tmp ./rmats/CONTROL_vs_${name%_rmats.csv}_tmp
done

conda deactivate

###END RMATS###


###START DESEQ2 and PLOT CREATION###
mkdir ./nextflow_results/DESEQ2

#Need to make DESEQ2 metadata sheet
awk -f $reformat_metadata_file samplesheet.csv > ./nextflow_results/DESEQ2/samplesheet_DESeq2.csv

conda activate R_nextflow_rnaseq

#Execute R script with reference files in positional arguments from source variables
Rscript $DESEQ2 $tx2gene $gtf_gene_annotation_table

conda deactivate

###END DESEQ2 and PLOT CREATION###

#Test if run is completed successfully by checking log file and if so delete the work directory
 if grep -q "Pipeline completed successfully" nextflow_RNAseq.log; then   
 rm -r ./work; 
 fi

