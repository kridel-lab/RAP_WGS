#----------------------------------------------------------------------

#Merged_001_evaluating_overlap_between_callers.R
#@Karin Isaev

#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS")

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
  "ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#obtain overlap between somatic variants called by Mutect2 and Strelka

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#these files below are summaries of VCF files produced by each tool
#to which additional soft filters were applied
#note these VCF files were normalized after each tool was run

strelka = list.files("/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR/strelka_filtered", pattern="filtered_strelka_calls")
strelka = strelka[-which(str_detect(strelka, "indels"))]
mutect2 = list.files("/cluster/projects/kridelgroup/RAP_ANALYSIS/MUTECT2_selected_VCFs", pattern="filtered_mutect2_calls.bed")

summary = as.data.frame(matrix(ncol=4, nrow=27)) ; colnames(summary) = c("patient", "unique_mutect2", "unique_strelka", "unique_common")

for(i in 1:length(strelka)){
  s_f = strelka[i]
  m_f = mutect2[i]
  #make sure same patient
  s_f_pat = paste(unlist(strsplit(s_f, "_"))[1:6], collapse="_")
  m_f_pat = paste(unlist(strsplit(m_f, "_"))[1:6], collapse="_")
  pat = paste(unlist(strsplit(m_f, "_"))[1:6], collapse="_")
  pat = unlist(strsplit(pat, "f"))[1]

  if(m_f_pat == s_f_pat){
    #read in each file
    s_f = fread(paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR/strelka_filtered", s_f, sep="/"))
    m_f = fread(paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/MUTECT2_selected_VCFs", m_f, sep="/"))
    both = merge(s_f, m_f, by = c("ChromKey", "CHROM", "POS", "mut_id", "REF", "ALT"))
    #save merged summary stats
    summary[i,] = c(pat, length(unique(m_f$mut_id)), length(unique(s_f$mut_id)), length(unique(both$mut_id)))
    #write merged summary mutation file
    both = both[,c("CHROM", "POS")]
    both$end = both$POS
    write.table(both, file=paste("merged_MUTECT2_STRELKA/", pat, "merged_mutations.bed", sep="_"), quote=F, col.names=F, row.names=F, sep="\t")
    print(summary)
  }
}

summary
summary$perc_overlap_strelka = as.numeric(summary$unique_common)/as.numeric(summary$unique_strelka)
summary$perc_overlap_mutect2 = as.numeric(summary$unique_common)/as.numeric(summary$unique_mutect2)

write.table(summary, file=paste("merged_MUTECT2_STRELKA/", "merged_mutations_data", "stats.txt", sep="_"), quote=F, row.names=F, sep="\t")

pats = as.data.frame(summary$patient)
write.table(pats, file="patient_ids.txt", quote=F, row.names=F, col.names=F, sep="\t")
