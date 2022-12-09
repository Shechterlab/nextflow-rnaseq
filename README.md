# nextflow-rnaseq-deseq2 overview

#Run nextflow-rnaseq pipeline followed by DESEQ2 with EnhancedVolcano Output on Einstein's HPC 

#The genome reference files are located at: /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/

#The companion scripts (also included here) can be found at: /gs/gsfs0/users/shechter-lab/data/NGS/scripts/maron

#Nextflow RNA-seq is already installed as a module on the HPC: https://nf-co.re/rnaseq

# 1. Clone the nextflow-rnaseq-deseq2-enhancedvolcano pipeline

git clone https://github.com/Shechterlab/nextflow-rnaseq-deseq2.git


# 2. Install Conda envrionments

conda env create -f R_nextflow_rnaseq.yml

conda env create -f zlib_nextflow_rnaseq.yml

# 3. Install R packages by activating the R conda envrionment, opening R, and then running the following commands (you shouldn't have to update the packages so can say no to that prompt)

conda activate R

R

#From within R run the below commands:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


BiocManager::install("DESeq2")

BiocManager::install("tximport")

BiocManager::install('EnhancedVolcano')

BiocManager::install("apeglm")

# 4. Once you have built the necessary envrionments, to run the pipeline first perform a nextflow test run 

sbatch ~/nextflow-rnaseq-deseq2/nextflow_rnaseq_test.sh

# 5. To run the pipeline, first make a directory for your experiment. Then, make a sample metadata file within that directory called: "nextflow_sample_metadata.csv" (example file here: https://github.com/nf-core/rnaseq/blob/3.9/assets/samplesheet.csv) with the experimental design. Example below:

mkdir ~/EB_2019_12_A549_test

nano ~/EB_2019_12_A549_test/nextflow_sample_metadata.csv

#past the following into the textfile, then save the changes

sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D4C1/D4C1_CRRA190011662-1a_HMKJGDSXX_L4_1.fq.gz,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D4C1/D4C1_CRRA190011662-1a_HMKJGDSXX_L4_2.fq.gz,unstranded
CONTROL_REP2,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D4C2/D4C2_CRRA190011663-1a_HMKJGDSXX_L4_1.fq.gz,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D4C2/D4C2_CRRA190011663-1a_HMKJGDSXX_L4_2.fq.gz,unstranded
CONTROL_REP3,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D4C3/D4C3_CRRA190011664-1a_HMKJGDSXX_L4_1.fq.gz,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D4C3/D4C3_CRRA190011664-1a_HMKJGDSXX_L4_2.fq.gz,unstranded
D2G_REP1,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D2G1/D2G1_CRRA190011644-1a_HMKJGDSXX_L2_1.fq.gz,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D2G1/D2G1_CRRA190011644-1a_HMKJGDSXX_L2_2.fq.gz,unstranded
D2G_REP2,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D2G2/D2G2_CRRA190011645-1a_HMKJGDSXX_L2_1.fq.gz,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D2G2/D2G2_CRRA190011645-1a_HMKJGDSXX_L2_2.fq.gz,unstranded
D2G_REP3,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D2G3/D2G3_CRRA190011646-1a_HMKJGDSXX_L4_1.fq.gz,/gs/gsfs0/users/shechter-lab/data/NGS/rna-seq/EB_2019_12_A549/fq/D2G3/D2G3_CRRA190011646-1a_HMKJGDSXX_L4_2.fq.gz,unstranded

# 6. Execute the pipeline

sbatch ~/nextflow-rnaseq-deseq2/nextflow_rnaseq_with_deseq2.sh ~/EB_2019_12_A549_test

#you can follow the output in the "~/nextflow-rnaseq-deseq2/nextflow_RNA_test.log" file or in the ~/EB_2019_12_A549_test/pipeline_info/ directory








