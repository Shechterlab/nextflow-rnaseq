# nextflow-rnaseq-deseq2 overview

#Run nextflow-rnaseq pipeline followed by DESEQ2 with EnhancedVolcano Output on Einstein's HPC 

#The genome reference files are located at: /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/

#The companion scripts (also included here) can be found at: /gs/gsfs0/users/shechter-lab/data/NGS/scripts/maron

#Nextflow RNA-seq is already installed as a module on the HPC: https://nf-co.re/rnaseq

# 1. Clone the nextflow-rnaseq-deseq2-enhancedvolcano pipeline

git clone https://github.com/Shechterlab/nextflow-rnaseq-deseq2.git


# 2. Install Conda envrionments

conda env create -f ~/nextflow-rnaseq-deseq2/R_nextflow_rnaseq.yml

conda env create -f ~/nextflow-rnaseq-deseq2/zlib_nextflow_rnaseq.yml

# 3. Install R packages by activating the R conda envrionment, opening R, and then running the following commands (you shouldn't have to update the packages so can say no to that prompt)

conda activate R_nextflow_rnaseq

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

# 5. To run the pipeline you only need a sample metadata file in your directory called "samplesheet.csv" (example file here: https://github.com/nf-core/rnaseq/blob/3.9/assets/samplesheet.csv) with the experimental design. Included example can be used as a test:

sbatch ~/nextflow-rnaseq-deseq2/nextflow_rnaseq_with_deseq2.sh

#you can follow the output in the "~/nextflow-rnaseq-deseq2/nextflow_RNAseq.log" file or in the ~/nextflow_results/pipeline_info/ directory








