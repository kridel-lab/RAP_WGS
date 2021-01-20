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
library(GenomicRanges)

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

strelka = list.files("/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR/strelka_filtered", pattern="indels.bed")
mutect2 = list.files("/cluster/projects/kridelgroup/RAP_ANALYSIS/MUTECT2_selected_VCFs", pattern="indels.bed")

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

    #get granges
    s_f$END=s_f$POS
    m_f$END=m_f$POS

    s_f_clean = s_f[,c("CHROM", "POS", "END", "ChromKey", "REF", "ALT")]
    s_f_clean$ChromKey = as.character(s_f_clean$ChromKey)
    s_f_clean$ChromKey = "*"
    colnames(s_f_clean)[2]="START"
    s_f_clean_gr = makeGRangesFromDataFrame(s_f_clean, keep.extra.columns=TRUE)

    m_f_clean = m_f[,c("CHROM", "POS", "END", "ChromKey", "REF", "ALT")]
    m_f_clean$ChromKey = as.character(m_f_clean$ChromKey)
    m_f_clean$ChromKey = "*"
    colnames(m_f_clean)[2]="START"
    m_f_clean_gr = makeGRangesFromDataFrame(m_f_clean, keep.extra.columns=TRUE)

    hits <- findOverlaps(s_f_clean_gr, m_f_clean_gr, ignore.strand=TRUE, maxgap=10)
    hits_overlap = cbind(s_f_clean[queryHits(hits),], m_f_clean[subjectHits(hits),])
    colnames(hits_overlap) = c("s_CHR", "s_START", "s_END", "s_STRAND",
    "s_REF", "s_ALT", "m_CHR", "m_START", "m_END", "m_STRAND", "m_REF", "s_ALT")
    hits_overlap$distance=hits_overlap$m_START-hits_overlap$s_START

    both = merge(s_f, m_f, by = c("ChromKey", "CHROM", "POS", "mut_id", "REF", "ALT"))
    #save merged summary stats
    summary[i,] = c(pat, length(unique(m_f$mut_id)), length(unique(s_f$mut_id)), length(unique(both$mut_id)))
    #write merged summary mutation file
    both = both[,c("CHROM", "POS")]
    both$end = both$POS
    write.table(both, file=paste("merged_MUTECT2_STRELKA/", pat, "merged_mutations_indels.bed", sep="_"), quote=F, col.names=F, row.names=F, sep="\t")
    print(summary)
    print(paste(dim(hits_overlap)[1], dim(both)[1]))
  }
}

summary
summary$perc_overlap_strelka = as.numeric(summary$unique_common)/as.numeric(summary$unique_strelka)
summary$perc_overlap_mutect2 = as.numeric(summary$unique_common)/as.numeric(summary$unique_mutect2)

write.table(summary, file=paste("merged_MUTECT2_STRELKA/", "merged_mutations_data", "stats_indels.txt", sep="_"), quote=F, row.names=F, sep="\t")

pats = as.data.frame(summary$patient)
write.table(pats, file="patient_ids.txt", quote=F, row.names=F, col.names=F, sep="\t")
