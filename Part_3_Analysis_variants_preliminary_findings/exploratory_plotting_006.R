#----------------------------------------------------------------------
#karin isaev
#last updated: June 24th 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("~/Documents/RAP_analysis")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats", "ggthemes")
lapply(packages, require, character.only = TRUE)

library(RColorBrewer)
library(openxlsx)
library(plotly)

display.brewer.all()
display.brewer.pal(9, "Set1")

#----------------------------------------------------------------------
#merging variants from Mutect2 and Strelka - summary 
#----------------------------------------------------------------------

res = fread("_merged_mutations_data_stats.txt")
res = res[order(-unique_mutect2, -unique_strelka, -unique_common)]
res$patient = factor(res$patient, levels=res$patient)
res = melt(res, id.vars = c("patient"), measure.vars = c("unique_mutect2", "unique_strelka", "unique_common"))
res$value = res$value/1000000
ggline(res, "patient", "value",
       linetype = "variable", shape = "variable",
       color = "variable", palette = brewer.pal(9, "Set1")[1:3]) + ylab("Mutations/Megabase") + grids(linetype = "dashed")+
  font("x.text", size = 8, color = "black") +
  rotate_x_text(55) + xlab("Sample ID")
ggsave(paste(date, "mutation_callers_comparison.pdf", sep="_"))

#as percentages 
res = fread("_merged_mutations_data_stats.txt")
res = res[order(-unique_mutect2, -unique_strelka, -unique_common)]
res$patient = factor(res$patient, levels=res$patient)
res = melt(res, id.vars = c("patient"), measure.vars = c("perc_overlap_strelka", "perc_overlap_mutect2"))
ggline(res, "patient", "value",
       linetype = "variable", shape = "variable",
       color = "variable", palette = brewer.pal(9, "Set1")[1:3]) + ylab("% of mutations called") + grids(linetype = "dashed")+
  font("x.text", size = 8, color = "black") +
  rotate_x_text(55) + xlab("Sample ID")
ggsave(paste(date, "mutation_callers_comparison.pdf", sep="_"))
ggsave(paste(date, "mutation_callers_comparison_percentages.pdf", sep="_"))

