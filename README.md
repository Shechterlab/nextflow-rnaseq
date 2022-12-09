# nextflow-rnaseq-deseq2 overview
Run nextflow-rnaseq pipeline followed by DESEQ2 with EnhancedVolcano Output on Einstein's HPC 
The genome reference files are located at: /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/
The companion scripts (also included here) can be found at: /gs/gsfs0/users/shechter-lab/data/NGS/scripts/maron
Nextflow RNA-seq is already installed as a module on the HPC: https://nf-co.re/rnaseq

# 1. Clone the nextflow-rnaseq-deseq2-enhancedvolcano pipeline

git clone https://github.com/Shechterlab/nextflow-rnaseq-deseq2.git


# 2. Install Conda envrionments

conda env create -f R.yml
conda env create -f zlib.yml

# 3. Install R packages by opening R and then running the following commands (you shouldn't have to update the packages so can say no to that prompt)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install('EnhancedVolcano')
BiocManager::install("apeglm")

# 4. To run the pipeline once it is installed, you need to make a sample metadata file called: "nextflow_sample_metadata.csv" (example file here: https://github.com/nf-core/rnaseq/blob/3.9/assets/samplesheet.csv)






