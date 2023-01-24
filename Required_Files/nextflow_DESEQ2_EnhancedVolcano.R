#This is intended to be run together with the nextflow rnaseq pipeline "nextflow_rnaseq.sh"
#Only required argument is the path to output directory for nextflow
#By Maxim Maron 12-28-22
#Load required packages
library(DESeq2)
library(tximport)
library(apeglm)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
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

#Set a directory where salmon files can be found
dir <- paste0(getwd(),'/nextflow_results/star_salmon')

#Create a variable to store the analysis name
analysis_name <- basename(getwd())

#Change to output directory for DESeq2
setwd(paste0(getwd(),'/nextflow_results/DESEQ2'))

#path to tx2gene file
tx2gene <- read.csv(args[1])

#Load in sample metadata
samples <- read.csv('samplesheet_DESeq2.csv', header = T)
#Load in quant files
files <- file.path(dir, samples$sample, "quant.sf")
names(files) <- samples$sample

#Load in Salmon quant files
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

#add new column to txi containing gene symbols derived from GTF using following command:
#zcat /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf.gz | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[5]"\t"$7}' | sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_biotype "//'| sed 's/gene_name "//' | sed 's/gene_biotype "//' | sed 's/"//g' | sed 's/ //g' | sed '1igene_id\tGeneSymbol\tChromosome\tClass\tStrand' > /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/hg38/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf_gene_annotation_table.txt

#path to gtf_gene_annotation_table.txt
features <- read.table(args[2], header =T)


#Build dds object
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ condition)

#filter out genes with less than 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#relevel the condition column so that the control is the reference
dds$condition <- relevel(dds$condition, ref = "CONTROL")

dds <- DESeq(dds)

# create a vst object from the dds object
vsd <- vst(dds, blind=FALSE)

# load the column data
cdata <- colData(dds)


# create a heatmap of the count matrix
pdf(paste0(analysis_name, "_Heatmap_CountMatrix.pdf"), height = 8, width = 8)
pheatmap(assay(vsd),
         cluster_rows = FALSE,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         annotation_col = as.data.frame(cdata[,"condition"], row.names=rownames(cdata)))
dev.off()

# calculate the distance between each sample
sampleDists <- dist(t(assay(vsd)))

# create a matrix of the sample distances
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste0(analysis_name, "_Heatmap_SampleDistances.pdf"), height = 8, width = 8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#Create a PCA of the sample data
pdf(paste0(analysis_name, "_PCA.pdf"), height = 8, width = 8)
plotPCA(vsd, intgroup=c("condition"))
dev.off()


# for loop to iterate through each comparison
for (i in 2:length(unique(resultsNames(dds)))){
  # name the comparison
  name <- unique(resultsNames(dds))[i]
  # run DESeq2 analysis
  res <- lfcShrink(dds, coef=name, type="apeglm")
  # add gene_id column to results
  res$gene_id <- rownames(res)
  # convert results to dataframe
  res_df <- as.data.frame(res)
  res_df <- merge(res_df, features, by.x="gene_id", by.y="gene_id")
  res_df <- res_df[order(res_df$padj),]
  write.csv(res_df, file=paste0(name,"_DESEQ2.csv"), row.names=F)
  pdf(paste0(name,"_MAplot.pdf"), height = 8, width = 8)
  plotMA(res)
  dev.off()
  res_df <- res_df[!is.na(res_df$padj),]
  pdf(paste0(name,'_Volcano_plot.pdf'), height = 8, width = 8)
  print(EnhancedVolcano(res_df,
                        lab = res_df$GeneSymbol,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = name,
                        pCutoff = 0.05,
                        FCcutoff = 1.5,
                        pointSize = 2.0,
                        labSize = 5.0,
                        col=c('#d4d4d4', '#d4d4d4', '#fdbb84', '#e34a33'),
                        colAlpha = 0.5))
  dev.off()
}

