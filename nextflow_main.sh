#!/bin/bash
#This script uses Nextflow to process RNA-seq data with nextflow then run DESEQ2 on the output
#This script requires a sample metadata file called "samplesheet.csv"
#e.g. Run as follows: "sbatch nextflow_rnaseq_with_deseq2.sh"

#By Maxim Maron on 12-09-22

#SBATCH -p unlimited      #partition/queue name
#SBATCH --job-name=nextflow_rnaseq    # Job name
#SBATCH --mail-type=NONE    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maxim.maron@einsteinmed.edu # Where to send mail
#SBATCH --ntasks=1      # Run on a single CPU
#SBATCH --cpus-per-task=40      # Number of CPU cores per task
#SBATCH --mem=80gb      # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=nextflow_RNAseq.log      # Standard output and error log
#SBATCH --get-user-env #Retrieve login envrionment variables

###RUNNING NEXTFLOW###

#Load conda envrionment and modules
.  /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate
conda activate zlib_nextflow_rnaseq
module load singularity
module load nextflow

#Export source variables to envrionment
.  /gs/gsfs0/users/shechter-lab/data/NGS/stds/source_files/20221220_nextflow_source_file.sh

nextflow run nf-core/rnaseq --input `pwd`/samplesheet.csv \
-profile singularity \
--outdir `pwd`/nextflow_results \
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
mkdir `pwd`/rmats/

#create an empty array called length that will store read lengths used to determine median read length for rmats
length=()

#go through all of the bam files and determine the read length, then append this value to the length array
#encountered error: "samtools: error while loading shared libraries: libncurses.so.5: cannot open shared object file: No such file or directory" that was fixed by going into the envrionment libraries and renaming libncurses.so.6 to libncurses.so.5

for file in $(ls `pwd`/nextflow_results/star_salmon/*.markdup.sorted.bam); do
	length+=( $(samtools view $file | head -n 1000 | gawk '{print length($10)}' | sort | uniq -c | perl -ane '$_ =~ s/^[ ]+//g;print $_' | sort -k 1nr,1nr | head -1 | cut -f2 -d " ") )
done

#determine the median read length and store this as a variable to be used with rmats
read_length=$(printf "%s\n" "${length[@]}" | datamash median 1)

#Make samplesheet that points to location of bam files for rmats
variables=$(awk -F',' 'FNR>1 && !a[$1]++{print $1}' `pwd`/samplesheet.csv)
for variable in $variables; do
	printf `pwd`/nextflow_results/star_salmon/${variable}.markdup.sorted.bam, >> `pwd`/rmats/${variable%%_REP*}_rmats_tmp.txt
	sed 's/\(.*\),/\1 /' `pwd`/rmats/${variable%%_REP*}_rmats_tmp.txt > `pwd`/rmats/${variable%%_REP*}_rmats.txt
done
rm `pwd`/rmats/*_rmats_tmp.txt


#build an array of file locations to analyze with rmats and remove the CONTROL
files_to_analyze=$(ls `pwd`/rmats/*.txt)

#run rmats using the array of variables
for file in $files_to_analyze; do
#only take the name after the last forward slash
	if [ $file != "`pwd`/rmats/CONTROL_rmats.txt" ]; then 
		name=$(echo $file | sed 's:.*/::')
		rmats.py --b1 `pwd`/rmats/CONTROL_rmats.txt --b2 $file --gtf $rmats_gtf -t paired --readLength $read_length --nthread 40 --od `pwd`/rmats/CONTROL_vs_${name%_rmats.txt} --tmp `pwd`/rmats/CONTROL_vs_${name%_rmats.txt}_tmp
	fi
done

#Remove temp directories
rm -r `pwd`/rmats/*_tmp
rm -r `pwd`/rmats/*/tmp
conda deactivate

###END RMATS###

###DESEQ2, CLUSTERPROFILER, AND PLOT CREATION###

#load conda R env
conda activate R_nextflow_rnaseq

#Identify directories with rmats analyses and then use Rscript to create violin plots and concatenated output files
directories=$(ls -d `pwd`/rmats/*/)
for directory in $directories; do
	Rscript $rmats_plots $directory
done

#Make directory for DESEQ2 analysis
mkdir `pwd`/nextflow_results/DESEQ2

#Need to make DESEQ2 metadata sheet
awk -f $reformat_metadata_file `pwd`/samplesheet.csv > `pwd`/nextflow_results/DESEQ2/samplesheet_DESeq2.csv

#Execute R script with reference files in positional arguments from source variables
Rscript $DESEQ2 $tx2gene $gtf_gene_annotation_table

#Make directory for CLUSTERPROFILER analysis
mkdir `pwd`/nextflow_results/ClusterProfiler

#Execute R script for Cluster Profiler
Rscript $ClusterProfiler

conda deactivate

###END DESEQ2 and PLOT CREATION###

#Test if run is completed successfully by checking log file and if so delete the work directory
 if grep -q "Pipeline completed successfully" `pwd`/nextflow_RNAseq.log; then   
 rm -r `pwd`/work; 
 fi
