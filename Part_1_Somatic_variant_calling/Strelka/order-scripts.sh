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
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_001_manta.sh

#----run strelka----------------------------------------------------------------

#call variants
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_002_strelka.sh
#select variants that pass filters
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_003_SelectVariants_post_strelka.sh
#annotate variants


#----run mutect2----------------------------------------------------------------

#1. first split normal samples by chromosomes
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_001_start_running_bam_splits_normal.sh

#2. split bam files by chromosomes
#first make sure there is enough space

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


#----Palimpsest-----------------------------------------------------------------

#prepare samples
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/prepare_annot_input_Palimpsest.R

#prepare CNAs
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/prepare_CNA_input_Palimpsest.R

#prepare SNVs
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/prepare_SNV_input_Palimpsest.R

#run palimpsest
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Palimpsest/Running_Palimpsest.R


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
