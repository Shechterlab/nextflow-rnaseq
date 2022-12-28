#Load libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(biomaRt)

#make a variable with the initial working directory path
initial_wd <- getwd()

analysis_name <- basename(getwd()) # define a variable with the base name of the working directory

setwd(paste0(getwd(),'/nextflow_results/DESEQ2')) # set the working directory to the nextflow_results/DESEQ2 directory

#Collect all files in my wd folder that end in ".MATS.JCEC.txt" and make a list of said files
file_list <- list.files(pattern = "_DESEQ2.csv$", recursive = TRUE)

#shorten name of each file to only reference the parent dataset
files <- gsub("_DESEQ2.csv", "", file_list)

#Read all files in "file_list"
dataset <- lapply(file_list, function(x){
  read.csv(x, header = TRUE)})

#Rename the dataframes in the list to match the original file from which they were derived
names(dataset) <- files

#Filter for padj < 0.05
dataset <- lapply(dataset, function(x) filter(x, padj < 0.05))
#Order by increasing padj
dataset <- lapply(dataset, function(x) {x <- x[order(x$padj,decreasing = F),];x})
#Filter for top 200 genes
dataset <- lapply(dataset, function(x) {x <- x[1:200,];x})

# load the dataset containing gene information from biomart
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://uswest.ensembl.org", dataset = 'hsapiens_gene_ensembl')

# get the entrez gene id for each gene
dataset_symbol <- lapply(dataset, function(x) getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_id"),values=x$gene_id,mart= mart))

# rename the column to entrez
dataset_symbol <- lapply(dataset_symbol, setNames, c("gene_id","entrez"))

# merge the entrez gene id with the original dataframe
dataframes <-mapply(merge, dataset, dataset_symbol, SIMPLIFY = FALSE)

# remove duplicates
dataframes <- lapply(dataframes, function(x) x[!duplicated(x$entrez),])

# order the dataframe by padj
dataframes <- lapply(dataframes, function(x) x[order(x$padj),])

# remove NA values
dataframes <- lapply(dataframes, function(x) x[!is.na(x$entrez),])

# bind the dataframes together
df <- dataframes %>% bind_rows(.id = 'df')

setwd(paste0(initial_wd,'/nextflow_results/ClusterProfiler')) # set the working directory to the nextflow_results/ClusterProfiler directory

# perform the GO enrichment analysis and output plots
formula_GO_BP <- compareCluster(entrez~df, data = df, fun="enrichGO",  OrgDb = org.Hs.eg.db, ont = 'BP', pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
formula_GO_simple_BP <-  simplify(formula_GO_BP, cutoff=0.7, by="p.adjust", select_fun=min)

pdf(paste0(analysis_name, "_enrichGO_BP.pdf"), height = 8, width = 8)
dotplot(formula_GO_simple_BP, showCategory = 20)+ scale_color_gradient(low = "#2c7bb6", high="#d7191c" ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

formula_GO_MF <- compareCluster(entrez~df, data = df, fun="enrichGO",  OrgDb = org.Hs.eg.db, ont = 'MF', pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
formula_GO_simple_MF <-  simplify(formula_GO_MF, cutoff=0.7, by="p.adjust", select_fun=min)

pdf(paste0(analysis_name, "_enrichGO_MF.pdf"), height = 8, width = 8)
dotplot(formula_GO_simple_MF, showCategory = 20)+ scale_color_gradient(low = "#2c7bb6", high="#d7191c" ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

formula_GO_CC <- compareCluster(entrez~df, data = df, fun="enrichGO",  OrgDb = org.Hs.eg.db, ont = 'CC', pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
formula_GO_simple_CC <-  simplify(formula_GO_CC, cutoff=0.7, by="p.adjust", select_fun=min)

pdf(paste0(analysis_name, "_enrichGO_CC.pdf"), height = 8, width = 8)
dotplot(formula_GO_simple_CC, showCategory = 20)+ scale_color_gradient(low = "#2c7bb6", high="#d7191c" ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()


formula_KEGG <- compareCluster(entrez~df, data = df, fun = 'enrichKEGG', pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

pdf(paste0(analysis_name, "_enrichKEGG.pdf"), height = 8, width = 8)
dotplot(formula_KEGG, showCategory = 20)+ scale_color_gradient(low = "#2c7bb6", high="#d7191c" ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

