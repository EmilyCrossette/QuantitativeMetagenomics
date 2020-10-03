###### SECTION 0: DOCUMENT SUMMARY  ---------------------------------------------------
# Marinobacter Recovery Analysis
# Emily Crossette
# 2020-08-27

# TO USE THIS SCRIPT: be sure to set working directory to folder containined data.

# This R script includes the computation of read recovery ratio across gene-length bins 
# to evaluate how target length biases read recovery using different read assignment 
# criteria. First, I mapped reads to a multi-fasta of individual genes and to the entire 
# spike-in DNA genome. I then wrote a combination of unix and python scripts to extract 
# gene-level read counts. This analysis takes those read counts and plots the read 
# recovery as a function of gene lengths. 

# This script was critical for developing our ultimate read-assignment approach to 
# limit biases caused by gene lengths and a precursor to python scripts we used to 
# generate a reproducible workflow. 

##### SECTION 1: LOAD PACKAGES, DATA & SETTINGS  -----------------------------------

library(tidyverse)
library(extrafont)
library(stringr)
library("binr")
library(Hmisc)
library(viridis)
library(scales)

# My favorite settings for ggplot
pretty_plot = theme_classic() + 
  theme(
    text = element_text(family = "Lucida Sans", color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    axis.line.y.left = element_line(color = "black"),
    #panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=15, face="bold"),
    axis.title=element_text(size=12,face="bold"),
    axis.text.y = element_text(size=12, color="#000000"),
    axis.text.x = element_text(size=12, color="#000000"),
    panel.spacing = unit(1.5, "lines")
    )

setwd("final-data")

Total_Reads = as.data.frame(
  read.table(
    file = "total_reads.csv", 
    sep =',', header = TRUE, stringsAsFactors = FALSE)) 
Total_Reads$Run_Sample = paste(Total_Reads$Run, Total_Reads$Sample, sep = "_")
Total_Reads = Total_Reads[c(-1,-2)]
Total_Reads$Read_Counts = as.numeric(Total_Reads$Read_Counts)

# Concentration in ng/uL of marinobacter calculated in spreadsheet 
genome_concentration = data.frame(
  Run_Sample = c("Run_1768_Sample_77304", "Run_1768_Sample_77305", "Run_1768_Sample_77306", "Run_1768_Sample_77307", "Run_1768_Sample_77308", "Run_1768_Sample_77309"),
  genome_conc = c(0.32876, 0.28386, 0.35719, 0.38334, 0.18098, 0.34862)
)

##### SECTION 2: Load Marinobactor Data from one Sample  -----------------------------------
# Reads mapping to multi-fasta of genes in paired kallisto setting
paired_gene = as.data.frame(
  read.table(file = "abund_GCA.tsv",
             sep ='\t', header = TRUE, quote = "", stringsAsFactors = FALSE))[c(1,4,2)]

# Reads mapping to genome in paired kallisto setting; extracted with bedtools
paired_genome = as.data.frame(
  read.table(file = "Sample_77304_CDS.txt",
             sep ='\t', header = FALSE, quote = "", stringsAsFactors = FALSE))[c(4,5,7)]

# Truncated average sequence depth 
paired_gene80 = as.data.frame(
  read.table(file = "TAD_CDS_77304.txt",
             sep ='\t', header = FALSE, quote = "", stringsAsFactors = FALSE))

# Reads mapping to multi-fasta of genes in unpaired paired kallisto setting
unpaired_genes = subset(
  as.data.frame(read.table(
    file = "MARINO_kalu.tsv",
    sep ='\t', header = TRUE, quote = "", stringsAsFactors = FALSE)),
  sid=="Sample_77304")[c(1,4,2)]

marino_length = unpaired_genes[c(1,3)]
  
#---------------------------------------------------------------#

##### SECTION 3: Merge Files Compute Gene Recovery Rate  -----------------------------------
#Merge all CARD and Marino Sets
colnames(paired_genome) = c("target_id", "est_counts", "length")
colnames(paired_gene80) = c("target_id", "est_counts", "length")

# The genes extracted from genome-mapping are counted as individual reads, not pairs. Need to multiply
# Gene-mapping counts by two to get number of reads (not read pairs)

paired_gene$est_counts = paired_gene$est_counts * 2

paired_genome$mapping = "genome-level"
paired_gene$mapping = "gene-level"
paired_gene80$mapping = "TAD80"
unpaired_genes$mapping = "un-paired, gene-level"

Sample_77304 = rbind(paired_gene, paired_genome)
rm(paired_gene, paired_genome)

Sample_77304 = rbind(Sample_77304, paired_gene80)
rm(paired_gene80)

Sample_77304 = rbind(Sample_77304, unpaired_genes)
rm(unpaired_genes)

# Calculate expected % abundance, using definition of Molecular Weight from: 
# https://www.idtdna.com/calc/Analyzer/Home/Definitions#MolecularWeight
# MW_A = 313.209, MW_C = 289.184, MW_G = 329.208, MW_T = 304.196

Avg_AT = (313.209 + 304.196)/2
Avg_GC = (289.184 + 329.208)/2

MW_genome_length = 4326849 
MW_genome_gc = 0.57
MW_genome = 2*((MW_genome_length * MW_genome_gc * Avg_GC) + 
                 (MW_genome_length - MW_genome_length*MW_genome_gc)*Avg_AT)

#Calculate gene copy concentration(copies/volume extract)
Read_recovery = Sample_77304 %>%
  dplyr::mutate(read_counts = 
                  if_else(mapping=="TAD80", est_counts*length/150, est_counts)) %>%
  dplyr::mutate(copies_uL =  # Copies of genes same as copies of genome
                  0.32876 * # ng/uL marinobactor *genome*
                  1E-9 * # g/ng
                  (1/MW_genome) * # moles/gram
                  6.022E23 # copies/mole
                  ) %>% 
  dplyr::mutate(Total_reads = 511522572,
                read_abund = read_counts / Total_reads, 
                read_abund_len = read_abund / length,
                observed_to_expected = read_abund_len/copies_uL)

#### SECTION 4: Binning genes by length and GC content for visualization #####

marino_length$length_quartile = with(marino_length, 
                                        cut(length, 
                                            breaks=quantile(length, 
                                                            probs=seq(0,1, by=0.05), 
                                                            na.rm=TRUE), 
                                            include.lowest=TRUE))


Marino_quantiles = marino_length[c(1,3)]
Read_recovery = merge(Read_recovery, Marino_quantiles, by = "target_id")

Marino_observed_len = Read_recovery %>%
  dplyr::select(mapping, length_quartile, length, observed_to_expected, target_id) %>%
  dplyr::group_by(length_quartile, mapping) %>%
  dplyr::summarise(min_len = min(length),
                  max_len = max(length),
                  mean_len = mean(length),
                  Quart_25  = quantile(observed_to_expected, 0.25),
                  Quart_75  = quantile(observed_to_expected, 0.75),
                  variability = IQR(observed_to_expected),
                  mean_obs_exp = mean(observed_to_expected),
                  num_genes = n_distinct(target_id),
                  standev = sd(observed_to_expected))

# plotcolor = c("#FE9F6DFF", "#DE4968FF", "#8C2981FF", "#3B0F70FF")
plotcolor = viridis_pal(option = "inferno", begin = 0.05, end = 0.7)(4)

# plotcolor = c("#4B2991", "#A3319F", "#EA4F88", "#FA7976")

# figure_1B = 
   ggplot(Marino_observed_len, aes(x= mean_len, y=mean_obs_exp, color=mapping, shape=mapping)) +
   geom_errorbar(aes(ymax=Quart_75, ymin=Quart_25)) +
   scale_color_manual(values=plotcolor) +
   scale_shape_manual(values=c(18, 15, 17, 16)) +
   geom_point(stat = "identity", size = 2) +
   ylim(0, 5E-14) +
   scale_x_continuous(trans='log2', breaks = c(250, 500, 1000, 2000, 4000), limits = c(240,4000)) +
   labs(
     x = "Average Gene Length of Bin",
     y = "Gene Recovery Ratio",
     fill = "Method") +
   pretty_plot +
   theme(axis.text.x = element_text(size=12, color="#000000", angle=-45, hjust=0.1),
         legend.title = element_text(size = 12),
         legend.text = element_text(size = 12),
         # legend.position = "none",
         # axis.title.y = element_blank(),
         # axis.text.y = element_blank()
         )
