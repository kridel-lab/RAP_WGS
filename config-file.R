#config file
#containing main folders and file names associated with scripts and programs
#@Karin isaev
#@July 20,2020

#folder containing mutect2 VCF files

#folder containing Strelka VCF files

#folder containing merged files from two callers
#that have been run through annovar

#folder containing VCF files and variant lists as input for Treeomics

#final matrix of SNVs merged with copy number data
#Our mutations
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

#CNAs annotated by gene that they overlap (no SNVs)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/TITAN_CNA/results/titan/hmm/optimalClusterSolution_files/titanCNA_ploidy2")
cnas_all = fread(fread(list.files(pattern="all_CNAs_protein_coding_samples.txt")[length(list.files(pattern="all_CNAs_protein_coding_samples.txt"))]))
