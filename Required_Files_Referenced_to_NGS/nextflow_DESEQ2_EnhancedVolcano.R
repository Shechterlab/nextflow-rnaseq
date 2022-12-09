#This is intended to be run together with the nextflow rnaseq pipeline "nextflow_rnaseq.sh"
#Only required argument is the path to output directory for nextflow
#By Maxim Maron 12-03-22
#Load required packages
library(DESeq2)
library(tximport)
library(apeglm)
library(EnhancedVolcano)
#Only load GenomicFeatures if making the tx2gene file for the first time
#library(GenomicFeatures)

#make script accept positiional command line arguments
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#Run code below once to make the tx2gene file for the appropriate gtf but then in subsequent steps just load in the csv because of memory limitations
#txdb <- makeTxDbFromGFF('/gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf.gz')
#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME")

setwd(paste0(args[1],'/DESEQ2'))

#path to tx2gene file
tx2gene <- read.csv(args[2])


dir <- paste0(args[1],'/star_salmon')
#Load in sample metadata
samples <- read.csv(paste0(args[1],'/DESEQ2/nextflow_sample_metadata.csv.DESEQ2.metadata.csv'), header = T)
#Load in quant files
files <- file.path(dir, samples$sample, "quant.sf")
names(files) <- samples$sample

#Load in Salmon quant files
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

#add new column to txi containing gene symbols derived from GTF using following command:
#zcat /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf.gz | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[5]"\t"$7}' | sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_biotype "//'| sed 's/gene_name "//' | sed 's/gene_biotype "//' | sed 's/"//g' | sed 's/ //g' | sed '1igene_id\tGeneSymbol\tChromosome\tClass\tStrand' > /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf_gene_annotation_table.txt

#path to gtf_gene_annotation_table.txt
features <- read.table(args[3], header =T)


#Build dds object
dds <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref = "CONTROL")

dds <- DESeq(dds)


for (i in 2:length(unique(resultsNames(dds)))){
  name <- unique(resultsNames(dds))[i]
  res <- lfcShrink(dds, coef=name, type="apeglm")
  res$gene_id <- rownames(res)
  res_df <- as.data.frame(res)
  res_df <- merge(res_df, features, by.x="gene_id", by.y="gene_id")
  res_df <- res_df[order(res_df$padj),]
  write.csv(res_df, file=paste0(name,"_DESEQ2.csv"), row.names=F)
  pdf(paste0(name,"_MAplot.pdf"), height = 8, width = 8)
  plotMA(res)
  dev.off()
  res_df <- res_df[!is.na(res_df$padj),]
  pdf(paste0(name,'_Vocano_plot.pdf'), height = 8, width = 8)
  print(EnhancedVolcano(res_df,
    lab = res_df$GeneSymbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    title = name,
    pCutoff = 0.05,
    FCcutoff = 1.5,
    pointSize = 2.0,
    labSize = 5.0,
    col=c('#d4d4d4', '#d4d4d4', '#fdbb84', '#e34a33'),
    colAlpha = 0.5))
  dev.off()
}
