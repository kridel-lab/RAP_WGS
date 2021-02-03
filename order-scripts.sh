#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_all_scripts

#----navigate to main directory-------------------------------------------------
cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#collect all sample BAM files into file
#ls */*/*recal.cram* > all_bam_files_raw.txt #with some manual clean up to only include tumour cram files

#collect all control samples
less all_control_samples.txt | wc -l

#----run manta------------------------------------------------------------------

#call SVs and indels
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_001_manta.sh

#----run strelka----------------------------------------------------------------

#1. call variants
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_002_strelka.sh

#2. select variants that pass filters
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_003_SelectVariants_post_strelka.sh

#3. select variants that pass filters (indels)
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_003_SelectVariants_post_strelka_indels.sh

#4. #mutationmapper
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/strelka_filtering_summary_mutationampper.R

#5. prepare bed files to merge with Mutect2
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/strelka_prepare_bed_files.R

#6. prepare bed files to merge with Mutect2 (indels)
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/strelka_prepare_bed_files_indels.R

#----run mutect2----------------------------------------------------------------

#1. first split normal samples by chromosomes
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_001_start_running_bam_splits_normal.sh

#2. split bam files by chromosomes
#requested temp dirctory for running mutect2
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_001_start_running_bam_splits.sh

#3. run mutect2 mutation calling on individual chr based bam files
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_001_running_mutect2.sh

#4. run getpileup on tumour samples
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_002A_running_GetPileupSummaries.sh

#5. run getpileup on normal samples
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_002A_running_GetPileupSummaries_control.sh

#6. get contamination tables
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_003A_running_CalculateContamination.sh

#7. run filter variant calls
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_004A_FilterVariants.sh

#8. run gather VCFs from across chromosomes
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_005A_running_gatherVCFs.sh

#9. select variants from combined VCF files to only include variants that passed
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_006_SelectVariants_post_mutect2.sh

#10. select variants from combined VCF files to only include variants that passed
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_006B_SelectVariants_post_mutect2_save_INDELS.sh

#11. soft filtering of mutect2 variants in R before merging with strelka
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_007_soft_filtering_mutect2_filtered_variants.R

#12. soft filtering of mutect2 variants in R before merging with strelka (indels)
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_007_soft_filtering_mutect2_filtered_variants_indels.R

#----merge variants from strelka and mutect2------------------------------------

#1. get list of mutations called by both tools
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_001_evaluating_overlap_between_callers.R

#2. get list of mutations called by both tools (indels)
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_001_evaluating_overlap_between_callers_indels.R

#3. use original mutect2 VCF files and keep only the mutations found by both callers
#run annovar
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_002_filter_variants_only_merged_run_annovar.sh

#4. use original mutect2 VCF files and keep only the mutations found by both callers (indels)
#run annovar
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_002_filter_variants_only_merged_run_annovar_indels.sh

#5. remove extra AF column from VCF files otherwise won't load properly into R
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_003B_remove_extra_AF_column.sh

#6. run annovar annotated VCF files through soft filtering
#vcfR functions don't work in 3.5.0
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_003_processing_merged_annovar.R

#7. run annovar annotated VCF files through soft filtering (indels)
#vcfR functions don't work in 3.5.0
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_003_processing_merged_annovar_indels.R

#8. some additional soft filtering and also the script that looks at overlaps
#between our gene mutations and several gene panels
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_004_additional_soft_filters_applied.R

#9. intersect mutations with CNAs and prepare for pyclone
#prepare mutations for bamreadcount
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_prepare_bed_files_to_fillter_VCFs_by_soft_filters_using_sequenza.R
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_prepare_bed_files_to_fillter_VCFs_by_soft_filters_indels_using_sequenza.R

#----TitanCNA-------------------------------------------------------------------

#run first run of TitanCNA
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/titan_cna_job.sh

#----Hatchet--------------------------------------------------------------------

#1. get BAM files for patients 001 and 002
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/Hatchet/001_p001_p002_get_BAM_files.sh

#2. run first part of Hatchet on patient 003 on the cluster (BAM files were already generated previously)
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/Hatchet/p003_hatchet_cluster_run.sh

#3. run first part of Hatchet on patient 001 on the cluster

#4. run first part of Hatchet on patient 002 on the cluster

#5. run second part of Hatchet on patient 001 locally using Gurobi

#6. run second part of Hatchet on patient 002 locally using Gurobi

#7. run second part of Hatchet on patient 003 locally using Gurobi

#8. prepare SNVs for each patient to be used for Hatchet provided explainMutationsCCF script
Rscript /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/Hatchet/prep_SNVs_for_CCF_calculation.R

#----CNVkit---------------------------------------------------------------------

#save control germline files as BAM files for all three samples in temp folder
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/CNVkit/CNVkit_pre_001_make_bam_files_for_control_samples.sh

#run WGS CNVkit mode on each tumour sample with appropriate control sample
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/CNVkit/CNVkit_001.sh

#run PureCN to get ploidy and purity estimates
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/CNVkit/CNVkit_002_PureCN.sh

#----Pyclone--------------------------------------------------------------------

#run bamreadcount
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Pyclone_get_missing_mutations.sh

#assemble results from bamreadcount into pyclone mutation list
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_008_get_pyclone_missing_counts.sh

#prepare input for pyclone-vi
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone-vi_001_make_input_files.R

#prepare input for pyclone-vi copy neutral assumption
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone-vi_001_make_input_files_neut.R

#prepare individual sample mutations
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_001_make_input_files.R

#run pyclone
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_002_run_main_program_all_muts.sh

#run pyclone-vi all muts
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone-vi_002_run_main_program_all_muts.sh

#run pyclone-vi all muts neutral assumption
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone-vi_002_run_main_program_all_muts_neut.sh

#run subset of mutations
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_002_run_main_program_subset_muts.sh

#run pyclone low mem
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_002_run_main_program_low_mem.sh

#prepare input for CITUP
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_003_make_citup_input.R

#run through CITUP
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_004_setting_up_CITUP.sh

#get output from CITUP
python /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_005_reading_CITUP_results_rap.py

#summarize results from pyclone
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_006_pyclone_summary.R

#run mapscape locally after results from citup have been transferred over
Rscript /Users/kisaev/github/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_007_mapscape.R

#----SciClone------------------------------------------------------------------

#generate input files
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/SciClone/sciclone_make_CNA_input_list.R
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/SciClone/sciclone_make_SNV_input_list.R

#run on individual samples
index=1 #from 1 to 20....
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/SciClone/sciclone_001_run.R $index

#run on all samples
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/SciClone/sciclone_002_job_submit.sh

#run pairwise sciclone comparing each sample to another sample 2D
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/SciClone/sciclone_003_job_submit.sh

#----Treeomics------------------------------------------------------------------

#1. make patient specific bed files from which to pull mutations in vcf files for treeomics
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_treeomics_prepare_bed_files_to_fillter_VCFs_by_soft_filters.R

#2. clean up VCF first for each sample to only inlcude mutations that want to include in analysis
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Treeomics/treeomics_001_A_input_VCF_files.sh

#3. change VCF file names in treeomics input to sample origin / organ
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Treeomics/treeomics_001_B_rename_sample_names_in_VCF_file.sh

#4. test run on small subset of protein coding genes
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Treeomics/treeomics_002_PCG_mutations_only_all.sh

#----Palimpsest-----------------------------------------------------------------

#prepare samples
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/prepare_annot_input_Palimpsest.R

#prepare CNAs
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/prepare_CNA_input_Palimpsest.R

#prepare SNVs
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/prepare_SNV_input_Palimpsest.R

#run palimpsest
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/Running_Palimpsest.R LY_RAP_0001
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/Running_Palimpsest.R LY_RAP_0002
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/Running_Palimpsest.R LY_RAP_0003

#----Summary-Plots--------------------------------------------------------------

#zoom in on DLBCL driver genes and how they are mutated
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_007_summary_driver_genes_across_samples.R


#----Staudt-classifier----------------------------------------------------------

#1. set up sample file
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Staudt_subtyping_samples/001_setting_up_sample_annotation_file.R

#2. set up mutation file
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Staudt_subtyping_samples/002_setting_up_mutation_file.R

#3. set up CNA file
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Staudt_subtyping_samples/003_setting_up_CNA_file.R

#4. summarize results from Staudt classifier
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Staudt_subtyping_samples/004_summarize_LymphGen_results.R

#DONE
