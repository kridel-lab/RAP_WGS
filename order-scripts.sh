#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_all_scripts

#----navigate to main directory-------------------------------------------------
/cluster/projects/kridelgroup/RAP_ANALYSIS

#----merge variants from strelka and mutect2------------------------------------

#get list of mutations called by both tools
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_001_evaluating_overlap_between_callers.R

#use original mutect2 VCF files and keep only the mutations found by both callers
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_002_filter_variants_only_merged_run_annovar.sh

#use original mutect2 VCF files and keep only the mutations found by both callers
#annotate this filtered VCF file using annovar
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_002_filter_variants_only_merged_run_annovar.sh

#remove extra AF column from VCF files otherwise won't load properly into R
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_003B_remove_extra_AF_column.sh

#run annovar annotated VCF files through soft filtering
#vcfR functions don't work in 3.5.0
module load R/3.6.1
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_003_processing_merged_annovar.R

#some additional soft filtering and also the script that looks at overlaps
#between our gene mutations and several gene panels
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_004_additional_soft_filters_applied.R

#----Pyclone--------------------------------------------------------------------

#prepare mutations for bamreadcount
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_prepare_bed_files_to_fillter_VCFs_by_soft_filters.R

#run bamreadcount
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Pyclone_get_missing_mutations.sh

#assemble results from bamreadcount into pyclone mutation list
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_008_get_pyclone_missing_counts.R

#prepare individual sample mutations
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_001_make_input_files.R

#run pyclone
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_002_run_main_program.sh

#run pyclone low mem
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_002_run_main_program_low_mem.sh

#prepare input for CITUP
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_003_make_citup_input.R

#run through CITUP
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_004_setting_up_CITUP.sh

#summarize results from pyclone


#----Treeomics------------------------------------------------------------------


#----PhyloWGS-------------------------------------------------------------------


#----Summary-Plots--------------------------------------------------------------
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_007_summary_driver_genes_across_samples.R
