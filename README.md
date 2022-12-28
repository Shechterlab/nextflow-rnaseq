# nextflow-rnaseq-deseq2 overview

#Run nextflow-rnaseq pipeline followed by rMATS and DESEQ2 with Splicing Event Violin Plots and EnhancedVolcano Output on Einstein's HPC 

#The genome reference files are located at: /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/

#The companion scripts (also included here) can be found at: /gs/gsfs0/users/shechter-lab/data/NGS/scripts/maron

#Nextflow RNA-seq is already installed as a module on the HPC: https://nf-co.re/rnaseq

# 1. Clone the nextflow-rnaseq-deseq2-enhancedvolcano pipeline 

#navigate to home directory and execute the code below in your home directory

cd

git clone https://github.com/Shechterlab/nextflow-rnaseq-deseq2.git


# 2. Install Conda envrionments

. ~/nextflow-rnaseq-deseq2/build_envs.sh

# 3. Install R packages 

#Activate the R conda envrionment, open R, and then running the following commands (you shouldn't have to update the packages so can say no to that prompt)

conda activate R_nextflow_rnaseq

R

#From within R run the below commands:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


BiocManager::install("DESeq2")

BiocManager::install("tximport")

BiocManager::install('EnhancedVolcano')

BiocManager::install("apeglm")

BiocManager::install("clusterProfiler")

BiocManager::install("org.Hs.eg.db")

BiocManager::install("biomaRt")

install.packages('tidyverse')

install.packages('data.table')

install.packages('RColorBrewer')

install.packages('pheatmap')


# 4. Perform a nextflow test run 

sbatch ~/nextflow-rnaseq-deseq2/nextflow_rnaseq_test.sh

# 5. Running the pipeline 

#You only need a sample metadata file in your directory called "samplesheet.csv" (example file here: https://github.com/nf-core/rnaseq/blob/3.9/assets/samplesheet.csv) with the experimental design. Please make sure the control samples within the samplesheet are clearly labeled "CONTROL_RepX". The included example can be used as a test and a guide.

sbatch ~/nextflow-rnaseq-deseq2/nextflow_rnaseq_with_deseq2.sh

#you can follow the output in the "nextflow_RNAseq.log" file or in the nextflow_results/pipeline_info/ directory








