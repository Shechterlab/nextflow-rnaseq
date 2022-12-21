#This is intended to be run together with the nextflow rnaseq pipeline "nextflow_main.sh"
#This script takes the rmats output from previous steps and outputs violin plots
#By Maxim Maron 12-21-22
#Load required packages
library(tidyverse)
library(data.table)

#make script accept positiional command line arguments
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#Set WD
setwd(args[1])
analysis_name <- basename(getwd())

#Collect all files in my wd folder that end in ".MATS.JCEC.txt" and make a list of said files
file_list <- list.files(pattern = ".MATS.JCEC.txt$", recursive = TRUE)

#shorten name of each file to only reference the parent dataset
files <- gsub(".MATS.JCEC.txt", "", file_list)

#Read all files in "file_list"
dataset <- lapply(file_list, function(x){
  read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)})

#Rename the dataframes in the list to match the original file from which they were derived
names(dataset) <- files

#Filter for FDR < 0.05
dataset <- lapply(dataset, function(x) filter(x, FDR < 0.05))

#Multiply by the RI dPSI by -1
dataset$RI$IncLevelDifference <- dataset$RI$IncLevelDifference * -1

#Combine the list into a single dataframe with a column indicating the parent dataset
df <- bind_rows(dataset, .id = 'df')

#Output this concatenated dataframe as a table
write.table(df, 'concatenated_rMATS_events_FDR_0.05.txt', quote = F, sep = '\t')

#Make a function to output summary statistics
stat_box_data <- function(y, upper_limit = 1) {
  return(
    data.frame(
      y = upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', median(y), '\n')
    )
  )
}


#Order the plot
order <- c("MXE","A3SS","A5SS","SE","RI")
df$df <- factor(df$df,  order)

#define colors for plot
colors<-c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e')

pdf(paste0(analysis_name,'_rMATS_violin_plots.pdf'), height = 11, width = 8)
ggplot(df, aes(x= IncLevelDifference,
               y=df,
               fill=df)) +
  geom_violin(alpha=0.5, orientation = 'y')+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values=colors)+
  #geom_text(aes(label = median), data = summary) +
  # scale_fill_manual(values=colors)+
  labs(title=paste0(analysis_name," ", "rMATS analysis \n"),
       x="deltaPSI",
       y = "Alternative Splicing Event",
       caption ="\n FDR <0.05 and RI deltaPSI * -1")+
  #geom_jitter(shape=21,width = .2, alpha = .2) +
  geom_vline(xintercept = 0 , linetype = 'dashed')+
  theme_classic()+
  theme(legend.position = 'na')+
  xlim(-1.1,1.1)+
  stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    hjust = 0.6,
    vjust = 1.2
  )
dev.off()
